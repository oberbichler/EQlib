#pragma once

#include "../Define.h"
#include "LinearSolverBase.h"

#include <Eigen/PardisoSupport>

namespace EQlib {

template <typename TSolver>
struct LinearSolverSetup
{
    static void apply(TSolver& solver)
    {
    }
};

template <>
struct LinearSolverSetup<Eigen::PardisoLDLT<Sparse, Eigen::Upper>>
{
    static void apply(Eigen::PardisoLDLT<Sparse, Eigen::Upper>& solver)
    {
        solver.pardisoParameterArray()[1] = 3;
    }
};

template <typename TSolver>
struct LinearSolver : LinearSolverBase
{
    TSolver m_solver;
    bool is_analyzed;

    LinearSolver() : is_analyzed(false)
    {
        LinearSolverSetup<TSolver>::apply(m_solver);
    }

    void analyze(Ref<const Sparse> a)
    {
        if (is_analyzed) {
            return;
        }

        m_solver.analyzePattern(a);

        is_analyzed = true;
    }

    void factorize(Ref<const Sparse> a)
    {
        analyze(a);
        m_solver.factorize(a);
    }

    void solve(Ref<const Vector> b, Ref<Vector> x)
    {
        x = m_solver.solve(b);
    }

    Eigen::ComputationInfo info() const
    {
        return m_solver.info();
    }
};

using SparseLU = LinearSolver<Eigen::SparseLU<Sparse>>;
using SymmetricConjugateGradient = LinearSolver<Eigen::ConjugateGradient<Sparse,
    Eigen::Upper>>;
using ConjugateGradient = LinearSolver<Eigen::ConjugateGradient<Sparse,
    Eigen::Lower | Eigen::Upper>>;
using PardisoLU = LinearSolver<Eigen::PardisoLU<Sparse>>;
using PardisoLLT = LinearSolver<Eigen::PardisoLLT<Sparse, Eigen::Upper>>;
using PardisoLDLT = LinearSolver<Eigen::PardisoLDLT<Sparse, Eigen::Upper>>;

} // namespace EQlib

