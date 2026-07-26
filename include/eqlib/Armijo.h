#pragma once

#include "Define.h"
#include "Log.h"
#include "Problem.h"
#include "Settings.h"
#include "Timer.h"

namespace eqlib {

class Armijo {
private: // types
    using Type = eqlib::Armijo;

private: // members
    Pointer<Problem> m_problem;
    double m_c;
    double m_rho;
    Vector m_x_init;
    Vector m_x;

public: // constructor
    Armijo(Pointer<Problem> problem)
        : m_problem(problem)
        , m_c(0.2)
        , m_rho(0.9)
        , m_x_init(problem->nb_variables())
        , m_x(problem->nb_variables())
    {
    }

public: // methods
    double search(Vector search_direction, double alpha_init, bool reset)
    {
        double alpha = alpha_init;

        m_x_init = m_problem->x();
        const double f_init = m_problem->f();
        const double cache = m_c * m_problem->df().dot(search_direction);

        m_x = m_x_init + alpha * search_direction;

        m_problem->set_x(m_x);
        m_problem->compute(0);
        double f = m_problem->f();

        while ((f - f_init) > (alpha * cache)) {
            alpha *= m_rho;
            m_x = m_x_init + alpha * search_direction;

            m_problem->set_x(m_x);
            m_problem->compute(0);
            f = m_problem->f();
        }

        if (reset) {
            m_problem->set_f(f_init);
            m_problem->set_x(m_x_init);
        }

        return alpha;
    }
}; // class Armijo

} // namespace eqlib
