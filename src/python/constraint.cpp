#include "common.h"

#include <eqlib/constraint.h>

void bind_constraint(py::module_ &m)
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

    bind<Constraint, PyConstraint>(m, "Constraint")
        .def(py::init<>())
        .def_property("is_active", &Constraint::is_active, &Constraint::set_is_active)
        .def_property("name", &Constraint::name, &Constraint::set_name)
        .def_property_readonly("nb_equations", &Constraint::nb_equations)
        .def_property_readonly("nb_parameters", &Constraint::nb_parameters)
        .def_property("equations", &Constraint::equations, &Constraint::set_equations)
        .def_property("parameters", &Constraint::parameters, &Constraint::set_parameters)
        .def_property("parameter_values", &Constraint::parameter_values, &Constraint::set_parameter_values);
}