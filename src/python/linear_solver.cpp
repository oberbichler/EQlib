#include "common.h"

#include <eqlib/linear_solver.h>

void bind_linear_solver(py::module &m, py::module &s)
{
    namespace py = pybind11;
    using namespace eqlib;
    using namespace pybind11::literals;

    struct PyLinearSolver : public LinearSolver {
        using LinearSolver::LinearSolver;

        virtual bool analyze(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a) override
        {
            pybind11::gil_scoped_acquire acquire;
            PYBIND11_OVERLOAD(bool, LinearSolver, analyze, ia, ja, a);
        }

        virtual bool factorize(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a) override
        {
            pybind11::gil_scoped_acquire acquire;
            PYBIND11_OVERLOAD(bool, LinearSolver, factorize, ia, ja, a);
        }

        virtual bool solve(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a, Ref<const Vector> b, Ref<Vector> x) override
        {
            pybind11::gil_scoped_acquire acquire;
            PYBIND11_OVERLOAD_PURE(bool, LinearSolver, factorize, ia, ja, a, b, x);
        }
    };

    py::class_<LinearSolver, PyLinearSolver, Pointer<LinearSolver>>(s, "LinearSolver")
        .def(py::init<>())
        .def_property("solver_name", &LinearSolver::solver_name, &LinearSolver::set_solver_name)
        .def("analyze", &LinearSolver::analyze, "ia"_a, "ja"_a, "a"_a)
        .def("factorize", &LinearSolver::factorize, "ia"_a, "ja"_a, "a"_a)
        .def("solve", &LinearSolver::solve, "ia"_a, "ja"_a, "a"_a, "b"_a, "x"_a);
}