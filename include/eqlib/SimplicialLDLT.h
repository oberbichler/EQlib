#pragma once

#include "Define.h"
#include "LinearSolver.h"

#include <Eigen/SparseCholesky>

#include <string>

namespace eqlib {

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

        Map<const ColMajorSparse> m(ia.size() - 1, ia.size() - 1, ja.size(), ia.data(), ja.data(), a.data());

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

        Map<const ColMajorSparse> m(ia.size() - 1, ia.size() - 1, ja.size(), ia.data(), ja.data(), a.data());

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

public: // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Base = LinearSolver;
        using Holder = Pointer<Type>;

        py::class_<Type, Base, Holder>(m, "SimplicialLDLT")
            // constructors
            .def(py::init<>());
    }
};

} // namespace eqlib