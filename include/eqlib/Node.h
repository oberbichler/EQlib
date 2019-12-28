#pragma once

#include "Define.h"
#include "Equation.h"
#include "Variable.h"

#include <string>

namespace eqlib {

class Node {
private: // types
    using Type = eqlib::Node;

private: // variables
    std::string m_name;

    Pointer<Variable> m_ref_x;
    Pointer<Variable> m_ref_y;
    Pointer<Variable> m_ref_z;
    Pointer<Variable> m_act_x;
    Pointer<Variable> m_act_y;
    Pointer<Variable> m_act_z;

    RobinMap<std::string, Pointer<Equation>> m_equations;
    RobinMap<std::string, Pointer<Variable>> m_variables;

public: // constructors
    Node(const double x, const double y, const double z) noexcept
        : m_ref_x(new_<Variable>(x))
        , m_ref_y(new_<Variable>(y))
        , m_ref_z(new_<Variable>(z))
        , m_act_x(new_<Variable>(x))
        , m_act_y(new_<Variable>(y))
        , m_act_z(new_<Variable>(z))
        , m_name("")
    {
    }

    Node() noexcept
        : Node(0, 0, 0)
    {
    }

public: // methods
    Pointer<Variable> ref_x() noexcept
    {
        return m_ref_x;
    }

    Pointer<Variable> ref_y() noexcept
    {
        return m_ref_y;
    }

    Pointer<Variable> ref_z() noexcept
    {
        return m_ref_z;
    }

    Pointer<Variable> x() noexcept
    {
        return m_act_x;
    }

    Pointer<Variable> y() noexcept
    {
        return m_act_y;
    }

    Pointer<Variable> z() noexcept
    {
        return m_act_z;
    }

    Vector3D ref_location() const noexcept
    {
        return Vector3D(m_ref_x->value(), m_ref_y->value(), m_ref_z->value());
    }

    void set_ref_location(Vector3D value) noexcept
    {
        m_ref_x->set_value(value[0]);
        m_ref_y->set_value(value[1]);
        m_ref_z->set_value(value[2]);
    }

    Vector3D act_location() const noexcept
    {
        return Vector3D(m_act_x->value(), m_act_y->value(), m_act_z->value());
    }

    void set_act_location(Vector3D value) noexcept
    {
        m_act_x->set_value(value[0]);
        m_act_y->set_value(value[1]);
        m_act_z->set_value(value[2]);
    }

    Vector3D displacements() const noexcept
    {
        return act_location() - ref_location();
    }

    void set_displacements(Vector3D value) noexcept
    {
        set_act_location(ref_location() + value);
    }

    const std::string& name() const
    {
        return m_name;
    }

    void set_name(const std::string& value)
    {
        m_name = value;
    }

    Pointer<Equation> equation(const std::string& name) noexcept
    {
        return m_equations[name];
    }

    bool has_equation(const std::string& name) noexcept
    {
        return m_equations.find(name) != m_equations.end();
    }

    Pointer<Variable> variable(const std::string& name) noexcept
    {
        if (name == "x") {
            return m_act_x;
        }
        if (name == "y") {
            return m_act_y;
        }
        if (name == "z") {
            return m_act_z;
        }
        if (name == "ref_x") {
            return m_ref_x;
        }
        if (name == "ref_y") {
            return m_ref_y;
        }
        if (name == "ref_z") {
            return m_ref_z;
        }

        if (m_variables.find(name) == m_variables.end()) {
            m_variables[name] = new_<Variable>();
        }

        return m_variables[name];
    }

    bool has_variable(const std::string& name) noexcept
    {
        return m_variables.find(name) != m_variables.end();
    }

public: // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
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
};

} // namespace eqlib