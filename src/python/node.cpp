#include "common.h"

#include <eqlib/node.h>

void bind_node(py::module_& m)
{
    namespace py = pybind11;
    using namespace eqlib;
    using namespace pybind11::literals;

    py::class_<Node, Pointer<Node>>(m, "Node", py::dynamic_attr())
        .def(py::init<>())
        .def(py::init<double, double, double, std::string>(), "x"_a = 0, "y"_a = 0, "z"_a = 0, "name"_a = "")
        .def_property("name", &Node::name, &Node::set_name)
        .def_property_readonly("x", &Node::x)
        .def_property_readonly("y", &Node::y)
        .def_property_readonly("z", &Node::z)
        .def_property("location", &Node::location, &Node::set_location);
}