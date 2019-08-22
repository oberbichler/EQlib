#pragma once

#include "../Define.h"
#include "LinearSolver.h"

#include <Eigen/PardisoSupport>

namespace EQlib {

template <typename TSolver>
struct EigenLinearSolver : LinearSolver
{
private:    // types
    template <typename TSolver>
    struct Setup
    {
        static void apply(TSolver& solver)
        {
        }
    };

    template <>
    struct Setup<Eigen::PardisoLDLT<Sparse, Eigen::Upper>>
    {
        static void apply(Eigen::PardisoLDLT<Sparse, Eigen::Upper>& solver)
        {
            solver.pardisoParameterArray()[1] = 3;
        }
    };

private:    // variables
    TSolver m_solver;
    bool is_analyzed;

public:     // constructors
    EigenLinearSolver() : is_analyzed(false)
    {
        Setup<TSolver>::apply(m_solver);
    }

public:     // methods
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

using SparseLU = EigenLinearSolver<Eigen::SparseLU<Sparse>>;
using SymmetricConjugateGradient = EigenLinearSolver<
    Eigen::ConjugateGradient<Sparse, Eigen::Upper>>;
using ConjugateGradient = EigenLinearSolver<Eigen::ConjugateGradient<Sparse,
    Eigen::Lower | Eigen::Upper>>;
using PardisoLU = EigenLinearSolver<Eigen::PardisoLU<Sparse>>;
using PardisoLLT = EigenLinearSolver<Eigen::PardisoLLT<Sparse, Eigen::Upper>>;
using PardisoLDLT = EigenLinearSolver<Eigen::PardisoLDLT<Sparse, Eigen::Upper>>;

} // namespace EQlib

