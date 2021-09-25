#include "common.h"

#include <eqlib/variable.h>

void bind_variable(py::module &m)
{
    using namespace eqlib;

    auto parameter_cls = bind<Variable>(m, "Variable");

    parameter_cls
        .def(py::init<const double, const bool, const std::string&>(), "value"_a = 0, "is_active"_a = true, "name"_a = "")
        .def_property("value", &Variable::value, &Variable::set_value)
        .def_property("is_active", &Variable::is_active, &Variable::set_is_active)
        .def_property("name", &Variable::name, &Variable::set_name)
        .def_property("lower_bound", &Variable::lower_bound, &Variable::set_lower_bound)
        .def_property("upper_bound", &Variable::upper_bound, &Variable::set_upper_bound);
}