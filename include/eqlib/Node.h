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
};

} // namespace eqlib
