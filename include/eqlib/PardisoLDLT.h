#pragma once

#include "Define.h"
#include "LinearSolver.h"

#include <mkl.h>

#include <array>
#include <string>
#include <vector>

namespace eqlib {

class PardisoLDLT : public LinearSolver {
private: // variables
    std::array<_MKL_DSS_HANDLE_t, 64> m_pt;
    MKL_INT m_mtype;
    std::array<MKL_INT, 64> m_iparm;
    std::vector<MKL_INT> m_perm;
    MKL_INT m_n;
    MKL_INT m_message_level;
    bool m_is_analyzed;

public: // constructors
    PardisoLDLT()
        : m_is_analyzed(false)
        , m_mtype(-2)
        , m_message_level(0)
    {
        ::pardisoinit(m_pt.data(), &m_mtype, m_iparm.data());

        m_iparm[0] = 1; // No solver default
        m_iparm[1] = 3; // 2: use Metis for the ordering, 3: parallel fill-in reducing ordering
        m_iparm[2] = 0; // Reserved. Set to zero. (??Numbers of processors, value of OMP_NUM_THREADS??)
        m_iparm[3] = 0; // No iterative-direct algorithm
        m_iparm[4] = 0; // No user fill-in reducing permutation
        m_iparm[5] = 0; // Write solution into x, b is left unchanged
        m_iparm[6] = 0; // Not in use
        m_iparm[7] = 2; // Max numbers of iterative refinement steps
        m_iparm[8] = 0; // Not in use
        m_iparm[9] = 13; // Perturb the pivot elements with 1E-13
        m_iparm[10] = 0; // 1: Use nonsymmetric permutation and scaling MPS
        m_iparm[11] = 0; // Not in use
        m_iparm[12] = 0; // 0: Maximum weighted matching algorithm is switched-off (default for symmetric).
        m_iparm[13] = 0; // Output: Number of perturbed pivots
        m_iparm[14] = 0; // Not in use
        m_iparm[15] = 0; // Not in use
        m_iparm[16] = 0; // Not in use
        m_iparm[17] = -1; // Output: Number of nonzeros in the factor LU
        m_iparm[18] = -1; // Output: Mflops for LU factorization
        m_iparm[19] = 0; // Output: Numbers of CG Iterations

        m_iparm[20] = 0; // 1x1 pivoting
        m_iparm[26] = 0; // No matrix checker
        m_iparm[27] = 0; // 0: double, 1: float
        m_iparm[34] = 1; // C indexing
        m_iparm[36] = 0; // CSR
        m_iparm[59] = 0; // 0 - In-Core ; 1 - Automatic switch between In-Core and Out-of-Core modes ; 2 - Out-of-Core
    }

private: // static methods
    static MKL_INT pardiso(
        _MKL_DSS_HANDLE_t pt,
        MKL_INT maxfct,
        MKL_INT mnum,
        MKL_INT mtype,
        MKL_INT phase,
        MKL_INT n,
        const double* a,
        const MKL_INT* ia,
        const MKL_INT* ja,
        MKL_INT* perm,
        MKL_INT nrhs,
        MKL_INT* iparm,
        MKL_INT msglvl,
        void* b,
        void* x)
    {
        MKL_INT error = 0;

        ::pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, b, x, &error);

        return error;
    }

public: // methods
    std::string solver_name() const override
    {
        return "MKL Pardiso";
    }

    bool analyze(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a) override
    {
        if (m_is_analyzed) {
            return false;
        }

        m_n = ia.size() - 1;
        m_perm.resize(m_n);

        MKL_INT error = pardiso(
            m_pt.data(), // pt
            1, // maxfct
            1, // mnum
            m_mtype, // mtype
            11, // phase
            m_n, // n
            a.data(), // a
            ia.data(), // ia
            ja.data(), // ja
            m_perm.data(), // perm
            0, // nrhs
            m_iparm.data(), // iparm
            m_message_level, // msglvl
            nullptr, // b
            nullptr // x
        );

        m_is_analyzed = true;

        return (error != 0);
    }

    bool factorize(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a) override
    {
        if (analyze(ia, ja, a)) {
            return true;
        }

        MKL_INT error = pardiso(
            m_pt.data(), // pt
            1, // maxfct
            1, // mnum
            m_mtype, // mtype
            22, // phase
            m_n, // n
            a.data(), // a
            ia.data(), // ia
            ja.data(), // ja
            m_perm.data(), // perm
            0, // nrhs
            m_iparm.data(), // iparm
            m_message_level, // msglvl
            nullptr, // b
            nullptr // x
        );

        return (error != 0);
    }

    bool solve(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a, Ref<const Vector> b, Ref<Vector> x) override
    {
        MKL_INT error = pardiso(
            const_cast<_MKL_DSS_HANDLE_t*>(m_pt.data()), // pt
            1, // maxfct
            1, // mnum
            m_mtype, // mtype
            33, // phase
            m_n, // n
            a.data(), // a
            ia.data(), // ia
            ja.data(), // ja
            const_cast<MKL_INT*>(m_perm.data()), // perm
            1, // nrhs
            const_cast<MKL_INT*>(m_iparm.data()), // iparm
            m_message_level, // msglvl
            const_cast<double*>(b.data()), // b
            x.data() // x
        );

        return (error != 0);
    }
};

} // namespace eqlib