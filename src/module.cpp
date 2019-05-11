#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <EQlib/Dof.h>
#include <EQlib/LBfgs.h>
#include <EQlib/Log.h>
#include <EQlib/LevenbergMarquardt.h>
#include <EQlib/NewtonDescent.h>
#include <EQlib/Node.h>
#include <EQlib/Parameter.h>
#include <EQlib/PyElement.h>
#include <EQlib/System.h>

template <typename Type, typename Module>
auto register_system(Module& m, std::string name)
{
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Holder = std::shared_ptr<Type>;

    return py::class_<Type, Holder>(m, name.c_str())
        // constructors
        .def(py::init<std::vector<std::shared_ptr<EQlib::Element>>, py::dict>(),
            "elements"_a, "linear_solver"_a = py::dict())
        // properties
        .def_property("load_factor", &Type::load_factor,
            &Type::set_load_factor)
        .def_property("x", &Type::x, &Type::set_x)
        // readonly properties
        .def_property_readonly("dofs", &Type::dofs)
        .def_property_readonly("nb_dofs", &Type::nb_dofs)
        .def_property_readonly("nb_elements", &Type::nb_elements)
        .def_property_readonly("elements", &Type::elements)
        .def_property_readonly("f", &Type::f)
        .def_property_readonly("g", &Type::g)
        .def_property_readonly("h", &Type::h)
        .def_property_readonly("message", &Type::message)
        .def_property_readonly("nb_free_dofs", &Type::nb_free_dofs)
        .def_property_readonly("nb_fixed_dofs", &Type::nb_fixed_dofs)
        .def_property_readonly("residual", &Type::residual)
        // methods
        .def("add_diagonal", &Type::add_diagonal, "value"_a)
        .def("compute", &Type::compute, "order"_a=2, "parallel"_a=true)
        .def("dof_index", &Type::dof_index)
        .def("element_indices", &Type::element_indices, "index"_a)
        .def("h_inv_v", &Type::h_inv_v)
        .def("h_v", &Type::h_v)
        .def("solve", &Type::solve, "maxiter"_a = 100, "rtol"_a = 1e-7,
            "xtol"_a = 1e-7, "regularization"_a=0.0, "parallel"_a=true)
        .def("solve_linear", &Type::solve_linear, "parallel"_a=true,
            "update_dofs"_a=true)
    ;
}

PYBIND11_MODULE(EQlib, m) {
    m.doc() = "EQlib by Thomas Oberbichler";
    m.attr("__author__") = "Thomas Oberbichler";
    m.attr("__copyright__") = "Copyright (c) 2018-2019, Thomas Oberbichler";
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
            .def("compute", &Type::compute, "g"_a, "h"_a)
        ;
    }

    // Log
    {
        using Type = EQlib::Log;

        py::class_<Type>(m, "Log")
            .def_property_static("info_level", [](py::object) { return
                Type::info_level(); }, [](py::object, const int value) {
                Type::set_info_level(value); })
            .def_static("debug", &Type::debug<const std::string&>, "message"_a)
            .def_static("info", py::overload_cast<const std::string&>(
                &Type::info<const std::string&>), "message"_a)
            .def_static("info", py::overload_cast<const int,
                const std::string&>(&Type::info<const std::string&>),
                "level"_a, "message"_a)
            .def_static("error", py::overload_cast<const std::string&>(
                &Type::error<const std::string&>), "message"_a)
            .def_static("error", py::overload_cast<const int,
                const std::string&>(&Type::error<const std::string&>),
                "level"_a, "message"_a)
            .def_static("warn", py::overload_cast<const std::string&>(
                &Type::warn<const std::string&>), "message"_a)
            .def_static("warn", py::overload_cast<const int,
                const std::string&>(&Type::warn<const std::string&>),
                "level"_a, "message"_a)
            .def_static("critical", py::overload_cast<const std::string&>(
                &Type::critical<const std::string&>), "message"_a)
            .def_static("critical", py::overload_cast<const int,
                const std::string&>(&Type::critical<const std::string&>),
                "level"_a, "message"_a)
        ;
    }

    // Node
    {
        using Type = EQlib::Node;
        using Holder = std::shared_ptr<Type>;

        py::class_<Type, Holder>(m, "Node", py::dynamic_attr())
            // constructors
            .def(py::init<>())
            .def(py::init<double, double, double>(), "x"_a=0, "y"_a=0, "z"_a=0)
            // readonly properties
            .def_property_readonly("x", &Type::x)
            .def_property_readonly("y", &Type::y)
            .def_property_readonly("z", &Type::z)
            // properties
            .def_property("ref_location", &Type::ref_location,
                &Type::set_ref_location)
            .def_property("act_location", &Type::act_location,
                &Type::set_act_location)
            .def_property("displacements", &Type::displacements,
                &Type::set_displacements)
            .def_property("forces", &Type::forces, &Type::set_forces)
            // methods
            .def("has_parameter", &Type::has_parameter, "name"_a)
            // operators
            .def("__getitem__", &Type::operator[],
                py::return_value_policy::reference_internal)
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
            .def_property("min_value", &Type::min_value, &Type::set_min_value)
            .def_property("max_value", &Type::max_value, &Type::set_max_value)
            .def_property("delta", &Type::delta, &Type::set_delta)
            .def_property("target", &Type::target, &Type::set_target)
            .def_property("result", &Type::result, &Type::set_result)
            .def_property("residual", &Type::residual, &Type::set_residual)
            .def_property("isfixed", &Type::isfixed, &Type::set_isfixed)
            .def_property("name", &Type::name, &Type::set_name)
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
            .def("__float__", [](const Type& self) { return self.act_value(); })
            .def("__repr__", &Type::to_string)
            .def("__copy__", [](const Type& self) { return Type(self); })
            .def("__deepcopy__", [](const Type& self, py::dict& memo) {
                return Type(self); }, "memodict"_a)
        ;
    }

    // System
    {
        using Type = EQlib::System<false>;

        register_system<Type>(m, "System");
    }

    // SymmetricSystem
    {
        using Type = EQlib::System<true>;

        register_system<Type>(m, "SymmetricSystem");
    }

    // Timer
    {
        using Type = EQlib::Timer;

        py::class_<Type>(m, "Timer")
            .def(py::init<>())
            .def("start", &Type::start)
            .def_property_readonly("ellapsed", &Type::ellapsed)
        ;
    }

    // LBfgs
    {
        using Type = EQlib::LBfgs;

        py::class_<Type>(m, "LBfgs")
            .def(py::init<std::shared_ptr<EQlib::System<true>>>(), "system"_a)
            .def("minimize", &Type::minimize, "maxiter"_a=100, "rtol"_a=1e-6,
                "xtol"_a=1e-6)
        ;
    }

    // LevenbergMarquardt
    {
        using Type = EQlib::LevenbergMarquardt;

        py::class_<Type>(m, "LevenbergMarquardt")
            .def(py::init<std::shared_ptr<EQlib::System<true>>>(), "system"_a)
            .def("minimize", &Type::minimize, "maxiter"_a=100, "rtol"_a=1e-6,
                "xtol"_a=1e-6)
        ;
    }

    // NewtonDescent
    {
        using Type = EQlib::NewtonDescent;

        py::class_<Type>(m, "NewtonDescent")
            .def(py::init<std::shared_ptr<EQlib::System<true>>>(), "system"_a)
            .def("minimize", &Type::minimize, "maxiter"_a=100, "rtol"_a=1e-6,
                "xtol"_a=1e-6, "line_search"_a = "none")
        ;
    }
}