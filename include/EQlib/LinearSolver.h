#pragma once

#include "Define.h"

#include <Eigen/PardisoSupport>

namespace EQlib {

class LinearSolver
{
private:    // variables
    Eigen::PardisoLDLT<Sparse, Eigen::Lower> m_solver;
    bool m_is_analyzed;

public:     // constructors
    LinearSolver() : m_is_analyzed(false)
    {
        m_solver.pardisoParameterArray()[1] = 3;    // enable OMP
    }

public:     // methods
    bool analyze(Ref<const Sparse> a)
    {
        if (m_is_analyzed) {
            return false;
        }

        m_solver.analyzePattern(a);

        m_is_analyzed = true;

        return (m_solver.info() != Eigen::Success);
    }

    bool factorize(Ref<const Sparse> a)
    {
        if (analyze(a)) {
            return true;
        }

        m_solver.factorize(a);

        return (m_solver.info() != Eigen::Success);
    }

    bool solve(Ref<const Vector> b, Ref<Vector> x) const
    {
        x = m_solver.solve(b);

        return (m_solver.info() != Eigen::Success);
    }
};

} // namespace EQlib