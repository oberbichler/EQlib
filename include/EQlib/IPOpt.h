#pragma once

#include "System.h"

#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>
#include <coin/IpTNLP.hpp>

#include <cassert>
#include <memory>

namespace EQlib {

class IPOptSystem : public Ipopt::TNLP
{
public:     // types
    using Index = Ipopt::Index;
    using Number = Ipopt::Number;

private:    // variables
    std::shared_ptr<System<true>> m_system;

public:     // constructors
    IPOptSystem(std::shared_ptr<System<true>> system)
    : m_system(system)
    { }

    ~IPOptSystem()
    { }

public:     // methods
    bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag,
        IndexStyleEnum& index_style)
    {
        // The problem described in IPOptSystem.hpp has n variables, x1, & x2,
        n = m_system->nb_free_dofs();

        // m equality constraint,
        m = 0;

        // 2 nonzeros in the jacobian (one for x1, and one for x2),
        nnz_jac_g = 0;

        // and 2 nonzeros in the hessian of the lagrangian
        // (one in the hessian of the objective for x2,
        //  and one in the hessian of the constraints for x1)
        nnz_h_lag = m_system->h().nonZeros();

        // We use the standard fortran index style for row/col entries
        index_style = C_STYLE;

        return true;
    }

    bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m,
        Number* g_l, Number* g_u)
    {
        // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
        // If desired, we could assert to make sure they are what we think they are.
        // assert(n == 2);
        // assert(m == 1);

        // x1 has a lower bound of -1 and an upper bound of 1
        for (int i = 0; i < m_system->nb_free_dofs(); i++) {
            const auto& dof = m_system->dof(i);
            x_l[i] = dof.lower_bound() - dof.ref_value();
            x_u[i] = dof.upper_bound() - dof.ref_value();
        }

        // x2 has no upper or lower bound, so we set them to
        // a large negative and a large positive number.
        // The value that is interpretted as -/+infinity can be
        // set in the options, but it defaults to -/+1e19
        // x_l[1] = -1.0e19;
        // x_u[1] = +1.0e19;

        // we have one equality constraint, so we set the bounds on this constraint
        // to be equal (and zero).
        // g_l[0] = g_u[0] = 0.0;

        return true;
    }

    bool get_starting_point(Index n, bool init_x, Number* x, bool init_z,
        Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda)
    {
        // Here, we assume we only have starting values for x, if you code
        // your own NLP, you can provide starting values for the others if
        // you wish.
        assert(init_x == true);
        assert(init_z == false);
        assert(init_lambda == false);

        // we initialize x in bounds, in the upper right quadrant
        for (int i = 0; i < m_system->nb_free_dofs(); i++) {
            const auto& dof = m_system->dof(i);
            x[i] = dof.delta();
        }

        return true;
    }

    bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
    {
        m_system->set_x(Map<const Vector>(x, n));
        m_system->assemble<0>(true);

        // return the value of the objective function

        obj_value = m_system->f();

        return true;
    }

    bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
    {
        m_system->set_x(Map<const Vector>(x, n));
        m_system->assemble<1>(true);

        // return the gradient of the objective function grad_{x} f(x)

        Map<Vector>(grad_f, n) = m_system->g();

        return true;
    }

    bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
    {
        // return the value of the constraints: g(x)
        // Number x1 = x[0];
        // Number x2 = x[1];

        // g[0] = -(x1 * x1 + x2 - 1.0);

        return true;
    }

    bool eval_jac_g(Index n, const Number* x, bool new_x, Index m,
        Index nele_jac, Index* iRow, Index* jCol, Number* values)
    {
        return true;
    }

    bool eval_h(Index n, const Number* x, bool new_x, Number obj_factor,
        Index m, const Number* lambda, bool new_lambda, Index nele_hess,
        Index* iRow, Index* jCol, Number* values
    )
    {
        if (values == nullptr) {
            const auto& h = m_system->h();

            int i = 0;

            for (int k = 0; k < h.outerSize(); ++k){
                for (Sparse::InnerIterator it(h, k); it; ++it){
                    iRow[i] = it.col();
                    jCol[i] = it.row();
                    i++;
                }
            }
        } else {
            m_system->set_x(Map<const Vector>(x, n));
            m_system->assemble<2>(true);

            const auto& h = m_system->h();

            int i = 0;

            for (int k = 0; k < h.outerSize(); ++k){
                for (Sparse::InnerIterator it(h, k); it; ++it){
                    values[i] = it.value();
                    i++;
                }
            }
        }

        return true;
    }

    void finalize_solution(Ipopt::SolverReturn status, Index n, const Number* x,
        const Number* z_L, const Number* z_U, Index m, const Number* g,
        const Number* lambda, Number obj_value, const Ipopt::IpoptData* ip_data,
        Ipopt::IpoptCalculatedQuantities* ip_cq)
    {
        // here is where we would store the solution to variables, or write to a file, etc
        // so we could use the solution. Since the solution is displayed to the console,
        // we currently do nothing here.
        for (int i = 0; i < m_system->nb_free_dofs(); i++) {
            const auto& dof = m_system->dof(i);
            dof.set_delta(x[i]);
        }
    }
};



class IPOpt
{
public:     // types
    using Index = Ipopt::Index;
    using Number = Ipopt::Number;

private:    // variables
    std::shared_ptr<System<true>> m_system;
    int m_info_level;
    int m_maxiter;
    double m_rtol;

public:     // constructors
    IPOpt(std::shared_ptr<System<true>> system)
    : m_info_level(0), m_system(std::move(system)), m_maxiter(100), m_rtol(1e-6)
    { }

public:     // methods
    int info_level() const
    {
        return m_info_level;
    }

    int maxiter() const
    {
        return m_maxiter;
    }

    void minimize()
    {
        // setup

        Log::info(1, "==> Minimizing nonlinear system...");
        Log::info(2, "Using IPOPT");

        Timer timer;

        Ipopt::SmartPtr<Ipopt::TNLP> problem = new IPOptSystem(m_system);

        Ipopt::SmartPtr<Ipopt::IpoptApplication> app =
            IpoptApplicationFactory();
        
        app->Options()->SetIntegerValue("max_iter", maxiter());
        app->Options()->SetIntegerValue("print_level", info_level());
        app->Options()->SetNumericValue("tol", rtol());

        Ipopt::ApplicationReturnStatus status = app->Initialize();

        if (status != Ipopt::Solve_Succeeded) {
            Log::error("Error during initialization");
        }

        status = app->OptimizeTNLP(problem);

        if (status == Ipopt::Solve_Succeeded) {
            // Retrieve some statistics about the solve
            Index iter_count = app->Statistics()->IterationCount();
            
            Log::info(2, "Problem solved in {} iterations", iter_count);

            Number final_obj = app->Statistics()->FinalObjective();

            Log::info(2, "Final value of the objective function is {}",
                final_obj);
        }

        Log::info(1, "System minimized in {:.3f} sec", timer.ellapsed());
    }

    double rtol() const
    {
        return m_rtol;
    }

    std::shared_ptr<System<true>> system() const
    {
        return m_system;
    }

    void set_info_level(int value)
    {
        m_info_level = value;
    }

    void set_maxiter(int value)
    {
        m_maxiter = value;
    }

    void set_rtol(double value)
    {
        m_rtol = value;
    }

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Type = IPOpt;

        py::class_<Type>(m, "IPOpt")
            // constructors
            .def(py::init<std::shared_ptr<EQlib::System<true>>>(), "system"_a)
            // read-only properties
            .def_property_readonly("system", &Type::system)
            // properties
            .def_property("info_level", &Type::info_level,
                &Type::set_info_level)
            .def_property("maxiter", &Type::maxiter,
                &Type::set_maxiter)
            .def_property("rtol", &Type::rtol, &Type::set_rtol)
            // methods
            .def("minimize", &Type::minimize)
        ;
    }
};

} // namespace EQlib