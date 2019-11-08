#pragma once

#include "Define.h"
#include "Problem.h"

#include <unsupported/Eigen/LevenbergMarquardt>

namespace EQlib {

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
        using JacobianType = Sparse;
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

        int df(const Vector& x, Sparse& fjac) const
        {
            m_problem->set_x(x);
            m_problem->compute<false, 1>();
            fjac = m_problem->dg();
            return 0;
        }

        index values() const
        {
            return m_problem->nb_dofs();
        }

        index inputs() const
        {
            return m_problem->nb_dofs();
        }
    };

public:     // constructor
    LevenbergMarquardt(const std::shared_ptr<Problem>& problem)
    : m_problem(problem)
    , m_iterations(0)
    , m_maxiter(0)
    , m_rnorm(0)
    , m_xnorm(0)
    , m_rtol(0)
    , m_xtol(0)
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

    void run(const index maxiter, const double rtol, const double xtol)
    {
        // setup

        Log::info(1, "==> Minimizing nonlinear system...");
        Log::info(2, "Using LM minimizer");

        Timer timer;

        Functor functor(m_problem);

        Vector x = m_problem->x();
        
        Eigen::LevenbergMarquardt<Functor> lm(functor);

        lm.setGtol(m_rtol);
        lm.setXtol(m_xtol);
        lm.setMaxfev(m_maxiter);
        
        lm.minimize(x);

        m_iterations = lm.iterations();
        m_rnorm = lm.gnorm();

        Log::info(1, "System minimized in {:.3f} sec", timer.ellapsed());
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

} // namespace EQlib