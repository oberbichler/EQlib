#include "common.h"

#include <eqlib/constraint.h>

void bind_constraint(py::module &m)
{
    namespace py = pybind11;
    using namespace eqlib;
    using namespace pybind11::literals;

    struct PyConstraint : public Constraint
    {
        using Constraint::Constraint;

        void compute(const Ref<Vector> g, const Ref<Matrix> dg, const Ref<Matrix> hm, const Request request) override
        {
            pybind11::gil_scoped_acquire acquire;
            PYBIND11_OVERLOAD_PURE(void, Constraint, compute, g, dg, hm, request);
        }
    };

    auto constraint_cls = bind<Constraint, PyConstraint>(m, "Constraint");

    constraint_cls
        .def(py::init<>())
        .def_property("is_active", &Constraint::is_active, &Constraint::set_is_active)
        .def_property("name", &Constraint::name, &Constraint::set_name)
        .def_property_readonly("nb_equations", &Constraint::nb_equations)
        .def_property_readonly("nb_variables", &Constraint::nb_variables)
        .def_property("equations", &Constraint::equations, &Constraint::set_equations)
        .def_property("variables", &Constraint::variables, &Constraint::set_variables)
        .def_property("variable_values", &Constraint::variable_values, &Constraint::set_variable_values);
}