#include "bindings.h"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <eqlib/Armijo.h>
#include <eqlib/NewtonRaphson.h>
#include <eqlib/SteepestDescent.h>
#include <eqlib/LinearSolver.h>
#include <eqlib/SimplicialLDLT.h>
#include <eqlib/SparseLU.h>
#ifdef EQLIB_USE_MKL
#include <eqlib/PardisoLDLT.h>
#endif

namespace eqlib {

void register_armijo(pybind11::module_& m)
{
    using Type = Armijo;
    namespace py = pybind11;
    using namespace pybind11::literals;

    py::class_<Type>(m, "Armijo")
        .def(py::init<Pointer<eqlib::Problem>>(), "problem"_a)
        // methods
        .def("search", &Type::search, py::call_guard<py::gil_scoped_release>(), "search_direction"_a, "alpha_init"_a = true, "reset"_a = true);
}

void register_newton_raphson(pybind11::module_& m)
{
    using Type = NewtonRaphson;
    namespace py = pybind11;
    using namespace pybind11::literals;

    py::class_<Type>(m, "NewtonRaphson")
        .def(py::init<Pointer<eqlib::Problem>>(), "problem"_a)
        .def("run", &Type::run, py::call_guard<py::gil_scoped_release>())
        // properties
        .def_property("damping", &Type::damping, &Type::set_damping)
        .def_property("maxiter", &Type::maxiter, &Type::set_maxiter)
        .def_property("rtol", &Type::rtol, &Type::set_rtol)
        // read-only properties
        .def_property_readonly("iterations", &Type::iterations)
        .def_property_readonly("rnorm", &Type::rnorm)
        .def_property_readonly("fevals", &Type::fevals)
        .def_property_readonly("gevals", &Type::gevals)
        .def_property_readonly("hevals", &Type::hevals);
}

void register_steepest_descent(pybind11::module_& m)
{
    using Type = SteepestDescent;
    namespace py = pybind11;
    using namespace pybind11::literals;

    py::class_<Type>(m, "SteepestDescent")
        .def(py::init<Pointer<eqlib::Problem>>(), "problem"_a)
        .def("run", &Type::run, py::call_guard<py::gil_scoped_release>())
        // properties
        .def_property("damping", &Type::damping, &Type::set_damping)
        .def_property("maxiter", &Type::maxiter, &Type::set_maxiter)
        .def_property("rtol", &Type::rtol, &Type::set_rtol)
        // read-only properties
        .def_property_readonly("iterations", &Type::iterations)
        .def_property_readonly("rnorm", &Type::rnorm)
        .def_property_readonly("fevals", &Type::fevals)
        .def_property_readonly("gevals", &Type::gevals)
        .def_property_readonly("hevals", &Type::hevals);

    // deprecated misspelled alias
    m.attr("SteepestDecent") = m.attr("SteepestDescent");
}

template <typename T>
class PyLinearSolver : public T {
public: // constructor
    using T::T;

public: // methods
    virtual bool analyze(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a) override
    {
        pybind11::gil_scoped_acquire acquire;
        PYBIND11_OVERLOAD(bool, T, analyze, ia, ja, a);
    }

    virtual bool factorize(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a) override
    {
        pybind11::gil_scoped_acquire acquire;
        PYBIND11_OVERLOAD(bool, T, factorize, ia, ja, a);
    }

    virtual bool solve(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a, Ref<const Vector> b, Ref<Vector> x) override
    {
        pybind11::gil_scoped_acquire acquire;
        PYBIND11_OVERRIDE_PURE(bool, T, factorize, ia, ja, a, b, x);
    }
};

void register_linear_solver(pybind11::module_& m)
{
    using Type = LinearSolver;
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Trampoline = PyLinearSolver<Type>;
    using Holder = Pointer<Type>;

    py::class_<Type, Trampoline, Holder>(m, "LinearSolver")
        .def(py::init<>())
        // read-only properties
        .def_property("solver_name", &Type::solver_name, &Type::set_solver_name)
        // methods
        .def("analyze", &Type::analyze, "ia"_a, "ja"_a, "a"_a)
        .def("factorize", &Type::factorize, "ia"_a, "ja"_a, "a"_a)
        .def("solve", &Type::solve, "ia"_a, "ja"_a, "a"_a, "b"_a, "x"_a);
}

void register_simplicial_ldlt(pybind11::module_& m)
{
    using Type = SimplicialLDLT;
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Base = LinearSolver;
    using Holder = Pointer<Type>;

    py::class_<Type, Base, Holder>(m, "SimplicialLDLT")
        // constructors
        .def(py::init<>());
}

void register_sparse_lu(pybind11::module_& m)
{
    using Type = SparseLU;
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Base = LinearSolver;
    using Holder = Pointer<Type>;

    py::class_<Type, Base, Holder>(m, "SparseLU")
        // constructors
        .def(py::init<>());
}

#ifdef EQLIB_USE_MKL

void register_pardiso_ldlt(pybind11::module_& m)
{
    using Type = PardisoLDLT;
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Base = LinearSolver;
    using Holder = Pointer<Type>;

    py::class_<Type, Base, Holder>(m, "PardisoLDLT")
        // constructors
        .def(py::init<>());
}

#endif // EQLIB_USE_MKL

} // namespace eqlib
