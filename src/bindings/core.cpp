#include "bindings.h"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <eqlib/Equation.h>
#include <eqlib/Variable.h>
#include <eqlib/Parameter.h>
#include <eqlib/Log.h>
#include <eqlib/Node.h>
#include <eqlib/Timer.h>
#include <eqlib/SparseStructure.h>

namespace eqlib {

void register_equation(pybind11::module_& m)
{
    using Type = Equation;
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Holder = Pointer<Type>;

    py::class_<Type, Holder>(m, "Equation")
        .def(py::init<double, double, double, std::string>(), "lower_bound"_a = -infinity, "upper_bound"_a = infinity, "multiplier"_a = 0.0, "name"_a = "")
        .def(py::init<>())
        .def_property("is_active", &Type::is_active, &Type::set_active)
        .def_property("lower_bound", &Type::lower_bound, &Type::set_lower_bound)
        .def_property("upper_bound", &Type::upper_bound, &Type::set_upper_bound)
        .def_property("name", &Type::name, &Type::set_name)
        .def_property("multiplier", &Type::multiplier, &Type::set_multiplier)
        .def("__repr__", &Type::to_string);
}

void register_variable(pybind11::module_& m)
{
    using Type = Variable;
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Holder = Pointer<Type>;

    py::class_<Type, Holder>(m, "Variable")
        // constructors
        .def(py::init<double, double, double, bool, double, std::string>(), "value"_a = 0.0, "lower_bound"_a = -infinity, "upper_bound"_a = infinity, "is_active"_a = true, "multiplier"_a = 1.0, "name"_a = "")
        .def(py::init<>())
        // methods
        .def("__float__", &Type::operator double)
        .def("clamp", &Type::clamp)
        // properties
        .def_property("value", py::overload_cast<>(&Type::value, py::const_), &Type::set_value)
        .def_property("lower_bound", &Type::lower_bound, &Type::set_lower_bound)
        .def_property("upper_bound", &Type::upper_bound, &Type::set_upper_bound)
        .def_property("is_active", &Type::is_active, &Type::set_active)
        .def_property("multiplier", &Type::multiplier, &Type::set_multiplier)
        .def_property("name", &Type::name, &Type::set_name);
}

void register_parameter(pybind11::module_& m)
{
    using Type = Parameter;
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Holder = Pointer<Type>;

    py::class_<Type, Holder>(m, "Parameter")
        .def(py::init<>())
        .def(py::init<double, std::string>(), "value"_a = 0.0, "name"_a = "")
        .def_property("value", &Type::value, &Type::set_value)
        .def_property("name", &Type::name, &Type::set_name)
        .def("__float__", &Type::operator double);
}

void register_log(pybind11::module_& m)
{
    using Type = Log;
    namespace py = pybind11;
    using namespace pybind11::literals;

    py::class_<Type>(m, "Log")
        .def_property_static("info_level", [](py::object) { return Type::info_level(); }, [](py::object, const int value) { Type::set_info_level(value); })
        .def_static("debug", [](const std::string& message) { Type::debug("{}", message); }, "message"_a)
        .def_static("info", [](const std::string& message) { Type::info("{}", message); }, "message"_a)
        .def_static("info", [](const int level, const std::string& message) { Type::info(level, "{}", message); },
            "level"_a, "message"_a)
        .def_static("error", [](const std::string& message) { Type::error("{}", message); }, "message"_a)
        .def_static("error", [](const int level, const std::string& message) { Type::error(level, "{}", message); },
            "level"_a, "message"_a)
        .def_static("warn", [](const int level, const std::string& message) { Type::warn(level, "{}", message); },
            "level"_a, "message"_a)
        .def_static("critical", [](const std::string& message) { Type::critical("{}", message); }, "message"_a)
        .def_static("critical", [](const int level, const std::string& message) { Type::critical(level, "{}", message); },
            "level"_a, "message"_a);
}

void register_node(pybind11::module_& m)
{
    using Type = Node;
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Holder = Pointer<Type>;

    py::class_<Type, Holder>(m, "Node", py::dynamic_attr())
        // constructors
        .def(py::init<>())
        .def(py::init<double, double, double>(), "x"_a = 0, "y"_a = 0, "z"_a = 0)
        // readonly properties
        .def_property_readonly("x", &Type::x)
        .def_property_readonly("y", &Type::y)
        .def_property_readonly("z", &Type::z)
        .def_property_readonly("ref_x", &Type::ref_x)
        .def_property_readonly("ref_y", &Type::ref_y)
        .def_property_readonly("ref_z", &Type::ref_z)
        // properties
        .def_property("ref_location", &Type::ref_location, &Type::set_ref_location)
        .def_property("act_location", &Type::act_location, &Type::set_act_location)
        .def_property("displacements", &Type::displacements, &Type::set_displacements)
        .def_property("name", &Type::name, &Type::set_name)
        // methods
        .def("equation", &Type::equation, "name"_a, py::return_value_policy::reference_internal)
        .def("variable", &Type::variable, "name"_a, py::return_value_policy::reference_internal)
        .def("has_equation", &Type::has_equation, "name"_a)
        .def("has_variable", &Type::has_variable, "name"_a);
}

void register_timer(pybind11::module_& m)
{
    using Type = Timer;
    namespace py = pybind11;
    using namespace py::literals;

    py::class_<Type>(m, "Timer")
        .def(py::init<>())
        .def("start", &Type::start)
        .def_property_readonly("elapsed", &Type::elapsed)
        // deprecated misspelled alias
        .def_property_readonly("ellapsed", &Type::elapsed);
}

template <typename TScalar, typename TIndex, bool TRowMajor, bool TIndexMap>
void register_sparse_structure(pybind11::module_& m, const std::string& name)
{
    using Type = SparseStructure<TScalar, TIndex, TRowMajor, TIndexMap>;
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Holder = Pointer<Type>;

    py::class_<Type, Holder>(m, name.c_str())
        // constructors
        .def(py::init<TIndex, TIndex, std::vector<TIndex>, std::vector<TIndex>>(), "rows"_a, "cols"_a, "ia"_a, "ja"_a)
        // static methods
        .def_static("convert_from", &Type::convert_from, "other"_a, "values"_a)
        .def_static("from_pattern", &Type::template from_pattern<std::vector<std::set<TIndex>>>, "rows"_a, "cols"_a, "pattern"_a)
        // methods
        .def("to_general", py::overload_cast<>(&Type::to_general, py::const_))
        .def("to_general", py::overload_cast<Ref<const Vector>>(&Type::to_general, py::const_))
        .def("get_index", &Type::get_index, "i"_a, "j"_a)
        .def("for_each", &Type::for_each, "action"_a)
        // read-only properties
        .def_property_readonly("rows", &Type::rows)
        .def_property_readonly("cols", &Type::cols)
        .def_property_readonly("nb_nonzeros", &Type::nb_nonzeros)
        .def_property_readonly("density", &Type::density)
        .def_property_readonly("ia", py::overload_cast<>(&Type::ia, py::const_))
        .def_property_readonly("ja", py::overload_cast<>(&Type::ja, py::const_));
}

void register_sparse_structures(pybind11::module_& m)
{
    register_sparse_structure<double, int, true, true>(m, "RowMajorSparseStructure");
    register_sparse_structure<double, int, false, true>(m, "ColMajorSparseStructure");
}

} // namespace eqlib
