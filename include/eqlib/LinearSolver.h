#pragma once

#include "Define.h"

#ifdef EQLIB_USE_MKL
#include <Eigen/PardisoSupport>
#else
#include <Eigen/SparseCholesky>
#endif

namespace eqlib {

class LinearSolver
{
private:    // variables
#ifdef EQLIB_USE_MKL
    Eigen::PardisoLDLT<Sparse, Eigen::Lower> m_solver;
#else
    Eigen::SimplicialLDLT<Sparse, Eigen::Lower> m_solver;
#endif
    bool m_is_analyzed;

public:     // constructors
    LinearSolver() : m_is_analyzed(false)
    {
#ifdef EQLIB_USE_MKL
        m_solver.pardisoParameterArray()[1] = 3;    // enable OMP
#endif
    }

public:     // static methods
    static std::string solver_name()
    {
#ifdef EQLIB_USE_MKL
        return "PardisoLDLT";
#else
        return "SimplicialLDLT";
#endif
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

} // namespace eqlib