#pragma once

#include "../Define.h"
#include "../Log.h"
#include "../Problem.h"
#include "../Settings.h"
#include "../Timer.h"

namespace EQlib {

class NewtonRaphson
{
private:    // types

private:    // members
    Pointer<Problem> m_problem;
    index m_iterations;
    index m_maxiter;
    index m_fevals;
    index m_gevals;
    index m_hevals;
    double m_rnorm;
    double m_xnorm;
    double m_rtol;
    double m_xtol;
    double m_damping;
    index m_stopping_reason;

private:    // methods

public:     // constructor
    NewtonRaphson(Pointer<Problem> problem)
    : m_problem(problem)
    , m_maxiter(100)
    , m_rtol(1e-6)
    , m_xtol(1e-6)
    , m_damping(0.0)
    {
        if (problem->is_constrained()) {
            throw std::runtime_error("Constraints are not supported");
        }
    }

public:     // methods
    index iterations() const noexcept
    {
        return m_iterations;
    }

    index fevals() const noexcept
    {
        return m_fevals;
    }

    index gevals() const noexcept
    {
        return m_gevals;
    }

    index hevals() const noexcept
    {
        return m_hevals;
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

    double damping() const noexcept
    {
        return m_damping;
    }

    void set_damping(const double value) noexcept
    {
        m_damping = value;
    }

    void run()
    {
        // setup

        Log::info(1, "==> Solving nonlinear system...");

        Timer timer;

        const index n = m_problem->nb_variables();

        for (index iteration = 0; ; iteration++) {
            // check max iterations

            if (iteration >= m_maxiter) {
                m_stopping_reason = 2;
                Log::info(2, "Stopped because iteration >= {}", m_maxiter);
                break;
            }

            Log::info(2, "Iteration {}", iteration + 1);

            // compute g and h

            Log::info(2, "Computing system...");

            m_problem->compute();
            m_gevals += 1;
            m_hevals += 1;

            Log::info(2, "The current value is {}", m_problem->f());

            // check residual

            Log::info(2, "Computing residual...");

            Vector m_residual = m_problem->df();

            const double rnorm = m_residual.norm();

            Log::info(2, "The norm of the residual is {}", rnorm);

            // check residual norm

            if (rnorm < m_rtol) {
                m_stopping_reason = 0;
                Log::info(2, "Stopped because rnorm < {}", m_rtol);
                break;
            }

            // solve iteration

            Log::info(2, "Solving the linear equation system...");

            if (m_damping != 0.0) {
                m_problem->hl_add_diagonal(m_damping);
            }

            Vector delta = m_problem->hl_inv_v(m_residual);

            // update system

            Log::info(2, "Updating system...");

            m_problem->set_x(m_problem->x() - delta);

            // check x norm

            const double xnorm = delta.norm();

            Log::info(2, "The norm of the step is {}", xnorm);

            if (xnorm < m_xtol) {
                Log::info(2, "Stopped because xnorm < {}", m_xtol);
                m_stopping_reason = 1;
                break;
            }
        }

        switch (m_stopping_reason) {
        case 2:
            Log::warn("The maximum number of iterations has been reached after {:.3f} sec", timer.ellapsed());
            break;
        case 3:
            Log::error("An unknown error has occurred after {:.3f} sec", timer.ellapsed());
            break;
        default:
            Log::info(1, "System solved in {:.3f} sec", timer.ellapsed());
            break;
        }
    }

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Type = EQlib::NewtonRaphson;

        py::class_<Type>(m, "NewtonRaphson")
            .def(py::init<Pointer<EQlib::Problem>>(), "problem"_a)
            .def("run", &Type::run)
            // properties
            .def_property("damping", &Type::damping, &Type::set_damping)
            .def_property("maxiter", &Type::maxiter, &Type::set_maxiter)
            .def_property("rtol", &Type::rtol, &Type::set_rtol)
            // read-only properties
            .def_property_readonly("iterations", &Type::iterations)
            .def_property_readonly("rnorm", &Type::rnorm)
            .def_property_readonly("fevals", &Type::fevals)
            .def_property_readonly("gevals", &Type::gevals)
            .def_property_readonly("hevals", &Type::hevals)
        ;
    }
};

} // namespace EQlib