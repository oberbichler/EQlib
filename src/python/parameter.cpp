#include "common.h"

#include <eqlib/parameter.h>

void bind_parameter(py::module_ &m)
{
    using namespace eqlib;

    bind<Parameter>(m, "Parameter")
        .def(py::init<const double, const bool, const std::string&>(), "value"_a = 0, "is_active"_a = true, "name"_a = "")
        .def(py::init<const double, std::pair<double, double>, const bool, const std::string&>(), "value"_a, "bounds"_a, "is_active"_a = true, "name"_a = "")
        .def_property("value", &Parameter::value, &Parameter::set_value)
        .def_property("is_active", &Parameter::is_active, &Parameter::set_is_active)
        .def_property("name", &Parameter::name, &Parameter::set_name)
        .def_property("lower_bound", &Parameter::lower_bound, &Parameter::set_lower_bound)
        .def_property("upper_bound", &Parameter::upper_bound, &Parameter::set_upper_bound);
}