#pragma once

#include "Define.h"
#include "LinearSolver.h"

#include <Eigen/SparseCholesky>

#include <string>

namespace eqlib {

class SimplicialLDLT : public LinearSolver
{
private:    // types
    using ColMajorSparse = Eigen::SparseMatrix<double, Eigen::ColMajor>;

private:    // variables
    Eigen::SimplicialLDLT<ColMajorSparse, Eigen::Lower> m_solver;
    bool m_is_analyzed;

public:     // constructors
    SimplicialLDLT() : m_is_analyzed(false)
    {
    }

public:     // methods
    std::string solver_name() const override
    {
        return "Eigen Simplicial LDLT";
    }

    bool analyze(Ref<const Sparse> a) override
    {
        if (m_is_analyzed) {
            return false;
        }

        Map<const ColMajorSparse> m(a.rows(), a.cols(), a.nonZeros(), a.outerIndexPtr(), a.innerIndexPtr(), a.valuePtr());

        m_solver.analyzePattern(m);

        const bool success = (m_solver.info() == Eigen::Success);

        if (success) {
            m_is_analyzed = true;
        }

        return !success;
    }

    bool factorize(Ref<const Sparse> a) override
    {
        if (analyze(a)) {
            return true;
        }

        Map<const ColMajorSparse> m(a.rows(), a.cols(), a.nonZeros(), a.outerIndexPtr(), a.innerIndexPtr(), a.valuePtr());

        m_solver.factorize(m);

        const bool success = (m_solver.info() == Eigen::Success);

        return !success;
    }

    bool solve(Ref<const Sparse> a, Ref<const Vector> b, Ref<Vector> x) override
    {
        x = m_solver.solve(b);

        const bool success = (m_solver.info() == Eigen::Success);

        return !success;
    }
};

} // namespace eqlib