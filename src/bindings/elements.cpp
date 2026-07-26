#include "bindings.h"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <eqlib/Constraint.h>
#include <eqlib/LambdaConstraint.h>
#include <eqlib/LambdaObjective.h>
#include <eqlib/Objective.h>

namespace eqlib {

template <typename T>
class PyConstraint : public T {
public: // constructor
    using T::T;

public: // methods
    void compute(Ref<Vector> rs, const std::vector<Ref<Vector>>& gs, const std::vector<Ref<Matrix>>& hs) const override
    {
        pybind11::gil_scoped_acquire acquire;
        PYBIND11_OVERRIDE_PURE(void, T, compute, rs, gs, hs);
    }
};

void register_constraint(pybind11::module_& m)
{
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Type = Constraint;
    using Trampoline = PyConstraint<Type>;
    using Holder = Pointer<Type>;

    py::class_<Type, Trampoline, Holder>(m, "Constraint")
        // constructors
        .def(py::init<>())
        .def(py::init<index, index>(), "nb_equations"_a, "nb_variables"_a)
        // read-only properties
        .def_property_readonly("nb_equations", &Type::nb_equations)
        .def_property_readonly("nb_variables", &Type::nb_variables)
        // properties
        .def_property("equations", &Type::equations, &Type::set_equations)
        .def_property("is_active", &Type::is_active, &Type::set_active)
        .def_property("name", &Type::name, &Type::set_name)
        .def_property("variables", &Type::variables, &Type::set_variables)
        // methods
        .def("compute", &Type::compute, "fs"_a, "gs"_a, "hs"_a)
        .def("equation", &Type::equation, "index"_a, py::return_value_policy::reference_internal)
        .def("variable", &Type::variable, "index"_a, py::return_value_policy::reference_internal);
}

template <typename T>
class PyObjective : public T {
public: // constructor
    using T::T;

public: // methods
    double compute(Ref<Vector> g, Ref<Matrix> h) const override
    {
        pybind11::gil_scoped_acquire acquire;
        PYBIND11_OVERRIDE_PURE(double, T, compute, g, h);
    }
};

void register_objective(pybind11::module_& m)
{
    using Type = Objective;
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Trampoline = PyObjective<Type>;
    using Holder = Pointer<Type>;

    py::class_<Type, Trampoline, Holder>(m, "Objective")
        // constructors
        .def(py::init<>())
        .def(py::init<index>(), "nb_variables"_a)
        // read-only properties
        .def_property_readonly("nb_variables", &Type::nb_variables)
        // properties
        .def_property("is_active", &Type::is_active, &Type::set_active)
        .def_property("name", &Type::name, &Type::set_name)
        .def_property("variables", &Type::variables, &Type::set_variables)
        // methods
        .def("compute", &Type::compute, "g"_a, "h"_a)
        .def("compute_all", [](const Type& self) {
            Vector g(self.nb_variables());
            Matrix h(self.nb_variables(), self.nb_variables());
            const double f = self.compute(g, h);
            return std::make_tuple(f, g, h);
        })
        .def("variable", &Type::variable, "index"_a, py::return_value_policy::reference_internal);
}

void register_lambda_constraint(pybind11::module_& m)
{
    using Type = LambdaConstraint;
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Base = Constraint;
    using Holder = Pointer<Type>;

    py::class_<Type, Base, Holder>(m, "LambdaConstraint")
        // constructors
        .def(py::init<const std::vector<Pointer<Equation>>&, const std::vector<Pointer<Variable>>&, Type::ComputeFunction>(), "equations"_a, "variables"_a, "compute"_a);
}

void register_lambda_objective(pybind11::module_& m)
{
    using Type = LambdaObjective;
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Base = Objective;
    using Holder = Pointer<Type>;

    py::class_<Type, Base, Holder>(m, "LambdaObjective")
        // constructors
        .def(py::init<const std::vector<Pointer<Variable>>&, Type::ComputeFunction>(), "variables"_a, "compute"_a);
}

} // namespace eqlib
