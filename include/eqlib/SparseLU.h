#pragma once

#include "Define.h"
#include "LinearSolver.h"
#include "SparseStructure.h"

#include <Eigen/SparseCholesky>

#include <string>

namespace eqlib {

class SparseLU : public LinearSolver {
private: // types
    using Type = SparseLU;
    using Sparse = Eigen::SparseMatrix<double, Eigen::ColMajor>;

private: // variables
    SparseStructure<double, int, true, true> m_structure;
    std::vector<index> m_value_indices;
    std::vector<double> m_a_values;
    Map<const Sparse> m_a;
    Eigen::SparseLU<Sparse> m_solver;
    bool m_is_analyzed;

public: // constructors
    SparseLU()
        : m_is_analyzed(false)
        , m_a(0, 0, 0, nullptr, nullptr, nullptr)
    {
        set_solver_name("EigenSparseLU");
    }

public: // methods
    bool analyze(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a) override
    {
        if (m_is_analyzed) {
            return false;
        }

        const int n = static_cast<int>(length(ia) - 1);

        const SparseStructure<double, int, true, true> symmetric_structure(n, n, ia, ja);

        const auto [structure, value_indices] = symmetric_structure.to_general();

        m_structure = structure;
        m_value_indices = value_indices;

        const auto nb_nonzeros = m_structure.nb_nonzeros();
        m_a_values.resize(nb_nonzeros);

        new (&m_a) Map<const Sparse>(n, n, nb_nonzeros, m_structure.ia().data(), m_structure.ja().data(), m_a_values.data());

        m_solver.analyzePattern(m_a);

        m_is_analyzed = true;

        return false;
    }

    bool factorize(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a) override
    {
        if (analyze(ia, ja, a)) {
            return true;
        }

        for (index i = 0; i < length(m_value_indices); i++) {
            m_a_values[i] = a[m_value_indices[i]];
        }

        m_solver.factorize(m_a);

        const bool success = (m_solver.info() == Eigen::Success);

        return !success;
    }

    bool solve(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a, Ref<const Vector> b, Ref<Vector> x) override
    {
        x = m_solver.solve(b.transpose());

        const bool success = (m_solver.info() == Eigen::Success);

        return !success;
    }

public: // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Base = LinearSolver;
        using Holder = Pointer<Type>;

        py::class_<Type, Base, Holder>(m, "SparseLU")
            // constructors
            .def(py::init<>());
    }
};

} // namespace eqlib