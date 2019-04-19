#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <EQlib/Dof.h>
#include <EQlib/Log.h>
#include <EQlib/Node.h>
#include <EQlib/Parameter.h>
#include <EQlib/PyElement.h>
#include <EQlib/System.h>

// #include <EQlib/IGA/LocationConstraint.h>
// #include <EQlib/IGA/Shell3D3P.h>

PYBIND11_MODULE(EQlib, m) {
    m.doc() = "EQlib by Thomas Oberbichler";
    m.attr("__author__") = "Thomas Oberbichler";
    m.attr("__copyright__") = "Copyright (c) 2019, Thomas Oberbichler";
    m.attr("__version__") = EQlib_VERSION;
    m.attr("__email__") = "thomas.oberbichler@gmail.com";
    m.attr("__status__") = "Development";

    namespace py = pybind11;
    using namespace pybind11::literals;

#if defined(EIGEN_USE_BLAS)
    m.attr("USE_BLAS") = true;
#else
    m.attr("USE_BLAS") = false;
#endif // EIGEN_USE_BLAS

#if defined(EIGEN_USE_MKL_ALL)
    m.attr("USE_MKL") = true;
#else
    m.attr("USE_MKL") = false;
#endif // EIGEN_USE_MKL_ALL

    // Dof
    {
        using Type = EQlib::Dof;

        py::class_<Type>(m, "Dof")
            .def_property("delta", &Type::delta, &Type::set_delta)
            .def_property("residual", &Type::residual, &Type::set_residual)
            .def_property_readonly("isfixed", &Type::isfixed)
            .def_property_readonly("target", &Type::target)
            .def("__eq__", &Type::operator==)
            .def("__hash__", &Type::hash)
        ;
    }

    // Element
    {
        using Type = EQlib::Element;
        using Trampoline = EQlib::PyElement;
        using Holder = std::shared_ptr<Type>;

        py::class_<Type, Trampoline, Holder>(m, "Element", py::dynamic_attr())
            .def(py::init<>())
            .def("dofs", &Type::dofs)
            .def("compute", &Type::compute)
        ;
    }

    // Log
    {
        using Type = EQlib::Log;

        py::class_<Type>(m, "Log")
            .def_property_static("info_level", [](py::object) { return
                Type::info_level; }, [](py::object, const int value) {
                Type::set_info_level(value); })
            .def_static("debug", &Type::debug<const std::string&>, "message"_a)
            .def_static("info", py::overload_cast<const std::string&>(
                &Type::info<const std::string&>), "message"_a)
            .def_static("info", py::overload_cast<const int,
                const std::string&>(&Type::info<const std::string&>),
                "level"_a, "message"_a)
            .def_static("error", &Type::error<const std::string&>, "message"_a)
            .def_static("warn", &Type::warn<const std::string&>, "message"_a)
            .def_static("critical", &Type::critical<const std::string&>,
                "message"_a)
        ;
    }

    // Node
    {
        using Type = EQlib::Node;
        using Holder = std::shared_ptr<Type>;

        py::class_<Type, Holder>(m, "Node", py::dynamic_attr())
            .def(py::init<>())
            .def(py::init<double, double, double>(), "x"_a=0, "y"_a=0, "z"_a=0)
            .def_property_readonly("x", &Type::x)
            .def_property_readonly("y", &Type::y)
            .def_property_readonly("z", &Type::z)
            .def_property("ref_location", &Type::ref_location,
                &Type::set_ref_location)
            .def_property("act_location", &Type::act_location,
                &Type::set_act_location)
            .def_property("displacements", &Type::displacements,
                &Type::set_displacements)
            .def_property("forces", &Type::forces, &Type::set_forces)
        ;
    }

    // Parameter
    {
        using Type = EQlib::Parameter;

        py::class_<Type>(m, "Parameter")
            .def(py::init<double, double, double, double, bool>(),
                "ref_value"_a, "act_value"_a, "target"_a=0, "result"_a=0,
                "isfixed"_a=false)
            .def(py::init<double, double, bool>(), "value"_a, "target"_a=0,
                "isfixed"_a=false)
            .def(py::init<>())
            .def_property("ref_value", &Type::ref_value, &Type::set_ref_value)
            .def_property("act_value", &Type::act_value, &Type::set_act_value)
            .def_property("delta", &Type::delta, &Type::set_delta)
            .def_property("target", &Type::target, &Type::set_target)
            .def_property("result", &Type::result, &Type::set_result)
            .def_property("residual", &Type::residual, &Type::set_residual)
            .def_property("isfixed", &Type::isfixed, &Type::set_isfixed)
            .def_property_readonly("dof", &Type::dof)
            .def(py::pickle([](const Type& self) {
                    return py::make_tuple(self.ref_value(), self.act_value(),
                        self.target(), self.result(), self.isfixed());
                }, [](py::tuple tuple) {
                    if (tuple.size() != 5) {
                        throw std::runtime_error("Invalid state!");
                    }

                    const auto ref_value = tuple[0].cast<double>();
                    const auto act_value = tuple[1].cast<double>();
                    const auto target = tuple[2].cast<double>();
                    const auto result = tuple[3].cast<double>();
                    const auto isfixed = tuple[4].cast<bool>();

                    return Type(ref_value, act_value, target, result, isfixed);
                }
            ))
            .def("__copy__", [](const Type& self) { return Type(self); })
            .def("__deepcopy__", [](const Type& self, py::dict& memo) {
                return Type(self); }, "memodict"_a)
        ;
    }

    // System
    {
        using Type = EQlib::System;

        py::class_<Type>(m, "System")
            .def(py::init<std::vector<std::shared_ptr<EQlib::Element>>>(),
                "elements"_a)
            .def_property_readonly("nb_dofs", &Type::nb_dofs)
            .def_property_readonly("nb_free_dofs", &Type::nb_free_dofs)
            .def_property_readonly("nb_fixed_dofs", &Type::nb_fixed_dofs)
            .def_property("nb_threads", &Type::nb_threads,
                &Type::set_nb_threads)
            .def_property_readonly("dofs", &Type::dofs)
            .def_property("load_factor", &Type::load_factor,
                &Type::set_load_factor)
            .def_property_readonly("lhs", &Type::lhs)
            .def_property_readonly("rhs", &Type::rhs)
            .def_property_readonly("residual", &Type::residual)
            .def_property_readonly("elements", &Type::elements)
            .def("compute", &Type::compute)
            .def("solve", &Type::solve)
            .def("solve_linear", &Type::solve_linear)
            .def_property_readonly("message",
                &Type::message)
            .def_property("x", &Type::x, &Type::set_x)
        ;
    }

    // // LocationConstraint
    // {
    //     using Type = EQlib::LocationConstraint;
    //     using Base = EQlib::Element;
    //     using Holder = std::shared_ptr<Type>;

    //     py::class_<Type, Base, Holder>(m, "LocationConstraint")
    //         .def(py::init<std::vector<std::shared_ptr<EQlib::Node>>,
    //             EQlib::Matrix, EQlib::Vector3D, double>(), "nodes"_a,
    //             "shape_functions"_a, "target"_a, "penalty"_a=1)
    //     ;
    // }

    // // Shell3D3P
    // {
    //     using Type = EQlib::Shell3D3P;
    //     using Base = EQlib::Element;
    //     using Holder = std::shared_ptr<Type>;

    //     py::class_<Type, Base, Holder>(m, "Shell3D3P")
    //         .def(py::init<std::vector<std::shared_ptr<EQlib::Node>>,
    //             EQlib::Matrix, double, double, double, double>(), "nodes"_a,
    //             "shape_functions"_a, "thickness"_a, "young_modulus"_a,
    //             "poisson_ratio"_a, "weight"_a)
    //     ;
    // }
}