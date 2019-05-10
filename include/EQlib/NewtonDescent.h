#pragma once

#include "System.h"

#include <memory>

namespace EQlib {

class NewtonDescent
{
private:    // members
    std::shared_ptr<System<true>> m_system;

private:    // methods
    double linesearch(Ref<const Vector> x, const Ref<Vector> search_dir,
        const double alpha_init = 1.0)
    {
        const double c = 0.2;
        const double rho = 0.9;
        double alpha = alpha_init;

        m_system->set_x(x + alpha * search_dir);
        m_system->assemble<0>(true);
        double f = m_system->f();

        m_system->set_x(x);
        m_system->assemble<1>(true);
        const double f_in = m_system->f();
        const Vector grad = m_system->g();

        const double cache = c * grad.dot(search_dir);

        while(f > f_in + alpha * cache) {
            alpha *= rho;

            m_system->set_x(x + alpha * search_dir);
            m_system->assemble<0>(true);
            f = m_system->f();
        }

        return alpha;
    }

public:     // constructor
    NewtonDescent(std::shared_ptr<System<true>> system)
    : m_system(std::move(system))
    {
    }

public:     // method
    void minimize(const int maxiter, const double rtol, const double xtol)
    {
        // setup

        Log::info(1, "==> Minimizing nonlinear system...");
        Log::info(2, "Using Newton Descent minimizer");

        Timer timer;

        const int nb_dofs = m_system->nb_free_dofs();

        Vector target(nb_dofs);

        for (int i = 0; i < nb_dofs; i++) {
            target[i] = m_system->dof(i).target();
        }

        if (m_system->load_factor() != 1.0) {
            target *= m_system->load_factor();
        }

        Vector residual(nb_dofs);
        Vector delta_x(nb_dofs);

        for (int iteration = 1; ; iteration++) {
            // check max iterations

            if (iteration >= maxiter) {
                Log::info(2, "Stopped because iteration >= {}", maxiter);
                Log::warn("The maximum number of iterations has been reached");
                break;
            }

            Log::info(2, "Iteration {}", iteration + 1);

            m_system->assemble<2>(true);

            // compute rnorm

            residual = target - m_system->g();

            const double rnorm = residual.norm();

            Log::info(2, "The norm of the residual is {}", rnorm);

            // check residual norm

            if (rnorm < rtol) {
                Log::info(2, "Stopped because rnorm < {}", rtol);
                break;
            }

            // linear solve

            m_system->add_diagonal(1e-5);

            delta_x = m_system->h_inv_v(residual);

            // linesearch

            const double rate = linesearch(m_system->x(), delta_x);

            Log::info(2, "Linesearch rate = {}", rate);

            delta_x *= rate;

            m_system->set_x(m_system->x() + delta_x);

            // compute xnorm

            const double xnorm = delta_x.norm();

            // check xnorm

            Log::info(2, "The norm of the step is {}", xnorm);

            if (xnorm < xtol) {
                Log::info(2, "Stopped because xnorm < {}", xtol);
                break;
            }
        }

        Log::info(1, "System minimized in {:.3f} sec", timer.ellapsed());
    }
};

} // namespace EQlib