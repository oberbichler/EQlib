#pragma once

#include "Define.h"
#include "Problem.h"

#include <unsupported/Eigen/LevenbergMarquardt>

namespace eqlib {

class LevenbergMarquardt
{
private:    // type
    using Type = LevenbergMarquardt;

private:    // members
    Pointer<Problem> m_problem;
    index m_iterations;
    index m_maxiter;
    double m_rnorm;
    double m_xnorm;
    double m_rtol;
    double m_xtol;

    struct Functor
    {
        using Scalar = double;
        using Index = index;
        using InputType = Vector;
        using ValueType = Vector;
        using JacobianType = Eigen::SparseMatrix<double, Eigen::ColMajor>;
        using QRSolver = Eigen::SparseQR<JacobianType, Eigen::COLAMDOrdering<int>>;

        enum {
            InputsAtCompileTime = Eigen::Dynamic,
            ValuesAtCompileTime = Eigen::Dynamic
        };

        std::shared_ptr<Problem> m_problem;

        Functor(const std::shared_ptr<Problem>& system) : m_problem(system)
        {
        }

        int operator()(const Vector& x, Vector& fvec) const
        {
            m_problem->set_x(x);
            m_problem->compute<false, 0>();
            fvec = m_problem->g();
            return 0;
        }

        int df(const Vector& x, JacobianType& fjac) const
        {
            m_problem->set_x(x);
            m_problem->compute<false, 1>();
            fjac = m_problem->dg().transpose();
            return 0;
        }

        index values() const
        {
            return m_problem->nb_equations();
        }

        index inputs() const
        {
            return m_problem->nb_variables();
        }
    };

public:     // constructor
    LevenbergMarquardt(const std::shared_ptr<Problem>& problem)
    : m_problem(problem)
    , m_rtol(1e-6)
    , m_xtol(1e-6)
    , m_maxiter(100)
    , m_rnorm(0)
    , m_xnorm(0)
    , m_iterations(0)
    {
    }

public:     // methods
    index iterations() const noexcept
    {
        return m_iterations;
    }

    index maxiter() const noexcept
    {
        return m_maxiter;
    }

    void set_maxiter(const index value) noexcept
    {
        m_maxiter = value;
    }

    double rnorm() const noexcept
    {
        return m_rnorm;
    }

    double rtol() const noexcept
    {
        return m_rtol;
    }

    void set_rtol(const double value) noexcept
    {
        m_rtol = value;
    }

    double xtol() const noexcept
    {
        return m_xtol;
    }

    void set_xtol(const double value) noexcept
    {
        m_xtol = value;
    }

    void run()
    {
        // setup

        m_rnorm = 0;
        m_xnorm = 0;
        m_iterations = 0;

        Log::task_begin("Solving nonlinear system using Levenberg-Marquardt...");

        Timer timer;

        Functor functor(m_problem);

        Vector x = m_problem->x();
        
        Eigen::LevenbergMarquardt<Functor> lm(functor);

        lm.setGtol(m_rtol);
        lm.setXtol(m_xtol);
        lm.setMaxfev(m_maxiter);
        
        lm.minimize(x);

        m_rnorm = lm.gnorm();
        m_iterations = lm.iterations();

        Log::task_end("System solved in {:.3f} sec", timer.ellapsed());
    }

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        py::class_<Type>(m, "LevenbergMarquardt")
            .def(py::init<Pointer<eqlib::Problem>>(), "problem"_a)
            .def("run", &Type::run)
            // properties
            .def_property("maxiter", &Type::maxiter, &Type::set_maxiter)
            .def_property("rtol", &Type::rtol, &Type::set_rtol)
            .def_property("xtol", &Type::xtol, &Type::set_xtol)
            // read-only properties
            .def_property_readonly("iterations", &Type::iterations)
            .def_property_readonly("rnorm", &Type::rnorm)
        ;
    }
};

} // namespace eqlib