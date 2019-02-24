#pragma once

#include "Define.h"
#include "Solver.h"

#include <Eigen/PardisoSupport>

namespace EQlib {

class SolverLDLT : public Solver {
private:    // types
    using SparseSolver = Eigen::PardisoLDLT<Sparse, Eigen::Upper>;

private:    // variables
    SparseSolver m_solver;

public:     // methods
    void analyze_pattern(Ref<const Sparse> a) override {
        m_solver.analyzePattern(a);
    }

    void set_matrix(Ref<const Sparse> a) override {
        m_solver.factorize(a);
    }

    void solve(Ref<const Vector> b, Ref<Vector> x) override {
        x = m_solver.solve(b);
    }
};

} // namespace EQlib