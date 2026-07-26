#include "bindings.h"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <eqlib/objectives/IgaMembrane3PAD.h>
#include <eqlib/objectives/IgaNormalDistanceAD.h>
#include <eqlib/objectives/IgaPointDistance.h>
#include <eqlib/objectives/IgaPointDistanceAD.h>
#include <eqlib/objectives/IgaPointLocation.h>
#include <eqlib/objectives/IgaRotationCouplingAD.h>
#include <eqlib/objectives/IgaShell3PAD.h>

namespace eqlib {

void register_iga_normal_distance_ad(pybind11::module_& m)
{
    using Type = IgaNormalDistanceAD;
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Holder = Pointer<Type>;
    using Base = Objective;

    py::class_<Type, Base, Holder>(m, "IgaNormalDistanceAD")
        .def(py::init<std::vector<Pointer<Node>>, std::vector<Pointer<Node>>>(), "nodes_a"_a, "nodes_b"_a)
        .def("add", &Type::add, "shape_functions"_a, "shape_functions_b"_a, "weight"_a);
}

void register_iga_membrane_3p_ad(pybind11::module_& m)
{
    using Type = IgaMembrane3PAD;
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Holder = Pointer<Type>;
    using Base = Objective;

    py::class_<Type, Base, Holder>(m, "IgaMembrane3PAD")
        .def(py::init<std::vector<Pointer<Node>>, double, double, double>(), "nodes"_a, "thickness"_a, "youngs_modulus"_a, "poissons_ratio"_a)
        .def("add", &Type::add, "shape_functions"_a, "weight"_a);
}

void register_iga_point_distance(pybind11::module_& m)
{
    using Type = IgaPointDistance;
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Holder = Pointer<Type>;
    using Base = Objective;

    py::class_<Type, Base, Holder>(m, "IgaPointDistance")
        .def(py::init<std::vector<Pointer<Node>>, std::vector<Pointer<Node>>>(), "nodes_a"_a, "nodes_b"_a)
        .def("add", &Type::add, "shape_functions"_a, "shape_functions_b"_a, "weight"_a);
}

void register_iga_point_distance_ad(pybind11::module_& m)
{
    using Type = IgaPointDistanceAD;
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Holder = Pointer<Type>;
    using Base = Objective;

    py::class_<Type, Base, Holder>(m, "IgaPointDistanceAD")
        .def(py::init<std::vector<Pointer<Node>>, std::vector<Pointer<Node>>>(), "nodes_a"_a, "nodes_b"_a)
        .def("add", &Type::add, "shape_functions"_a, "shape_functions_b"_a, "weight"_a);
}

void register_iga_point_location(pybind11::module_& m)
{
    using Type = IgaPointLocation;
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Holder = Pointer<Type>;
    using Base = Objective;

    py::class_<Type, Base, Holder>(m, "IgaPointLocation")
        .def(py::init<std::vector<Pointer<Node>>>(), "nodes"_a)
        .def("add", &Type::add, "shape_functions"_a, "target"_a, "weight"_a);
}

void register_iga_rotation_coupling_ad(pybind11::module_& m)
{
    using Type = IgaRotationCouplingAD;
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Holder = Pointer<Type>;
    using Base = Objective;

    py::class_<Type, Base, Holder>(m, "IgaRotationCouplingAD")
        .def(py::init<std::vector<Pointer<Node>>, std::vector<Pointer<Node>>>(), "nodes_a"_a, "nodes_b"_a)
        .def("add", &Type::add, "shape_functions"_a, "shape_functions_b"_a, "axis"_a, "weight"_a);
}

void register_iga_shell_3p_ad(pybind11::module_& m)
{
    using Type = IgaShell3PAD;
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Holder = Pointer<Type>;
    using Base = Objective;

    py::class_<Type, Base, Holder>(m, "IgaShell3PAD")
        .def(py::init<std::vector<Pointer<Node>>, double, double, double>(), "nodes"_a, "thickness"_a, "youngs_modulus"_a, "poissons_ratio"_a)
        .def("add", &Type::add, "shape_functions"_a, "weight"_a);
}

} // namespace eqlib
