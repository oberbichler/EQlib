#include "common.h"

#include <eqlib/equation.h>

void bind_equation(py::module_ &m)
{
    using namespace eqlib;

    auto equation_cls = bind<Equation>(m, "Equation");

    equation_cls
        .def(py::init<const bool, const std::string>(), "is_active"_a = true, "name"_a = "")
        .def_property("value", &Equation::value, &Equation::set_value)
        .def_property("is_active", &Equation::is_active, &Equation::set_is_active)
        .def_property("name", &Equation::name, &Equation::set_name)
        .def_property("lower_bound", &Equation::lower_bound, &Equation::set_lower_bound)
        .def_property("upper_bound", &Equation::upper_bound, &Equation::set_upper_bound);
}