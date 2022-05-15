#pragma once

#include "../common.h"
#include "../linear_solver.h"

#include <Eigen/SparseCholesky>

#include <string>

namespace eqlib {
namespace linear_solvers {

class SimplicialLDLT : public LinearSolver {
private: // types
    using Type = SimplicialLDLT;
    using ColMajorSparse = Eigen::SparseMatrix<double, Eigen::ColMajor>;

private: // variables
    Eigen::SimplicialLDLT<ColMajorSparse, Eigen::Lower> m_solver;
    bool m_is_analyzed;

public: // constructors
    SimplicialLDLT()
        : m_is_analyzed(false)
    {
        set_solver_name("SimplicialLDLT");
    }

public: // methods
    bool analyze(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a) override
    {
        if (m_is_analyzed) {
            return false;
        }

        const index size = ia.size() - 1;
        const index nnz = ja.size();

        Map<const ColMajorSparse> m(size, size, nnz, ia.data(), ja.data(), a.data());

        m_solver.analyzePattern(m);

        const bool success = (m_solver.info() == Eigen::Success);

        if (success) {
            m_is_analyzed = true;
        }

        return !success;
    }

    bool factorize(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a) override
    {
        if (analyze(ia, ja, a)) {
            return true;
        }

        const index size = ia.size() - 1;
        const index nnz = ja.size();

        Map<const ColMajorSparse> m(size, size, nnz, ia.data(), ja.data(), a.data());

        m_solver.factorize(m);

        const bool success = (m_solver.info() == Eigen::Success);

        return !success;
    }

    bool solve(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a, Ref<const Vector> b, Ref<Vector> x) override
    {
        x = m_solver.solve(b.transpose());

        const bool success = (m_solver.info() == Eigen::Success);

        return !success;
    }
};

} // namespace linear_solvers
} // namespace eqlib