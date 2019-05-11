#pragma once

#include "Define.h"
#include "Parameter.h"

#include <memory>
#include <string>
#include <unordered_map>

namespace EQlib {

class Node
{
private:    // variables
    Parameter m_x;
    Parameter m_y;
    Parameter m_z;

    std::unordered_map<std::string, Parameter> m_parameters;

public:     // constructors
    Node(const double x, const double y, const double z) noexcept
    : m_x(x)
    , m_y(y)
    , m_z(z)
    { }

    Node() noexcept : Node(0, 0, 0)
    { }

public:     // getters and setters
    Parameter& x() noexcept
    {
        return m_x;
    }

    Parameter& y() noexcept
    {
        return m_y;
    }

    Parameter& z() noexcept
    {
        return m_z;
    }

    Vector3D ref_location() const noexcept
    {
        return Vector3D(m_x.ref_value(), m_y.ref_value(), m_z.ref_value());
    }

    void set_ref_location(Vector3D value) noexcept
    {
        m_x.set_ref_value(value[0]);
        m_y.set_ref_value(value[1]);
        m_z.set_ref_value(value[2]);
    }

    Vector3D act_location() const noexcept
    {
        return Vector3D(m_x.act_value(), m_y.act_value(), m_z.act_value());
    }

    void set_act_location(Vector3D value) noexcept
    {
        m_x.set_act_value(value[0]);
        m_y.set_act_value(value[1]);
        m_z.set_act_value(value[2]);
    }

    Vector3D displacements() const noexcept
    {
        return Vector3D(m_x.delta(), m_y.delta(), m_z.delta());
    }

    void set_displacements(Vector3D value) noexcept
    {
        m_x.set_delta(value[0]);
        m_y.set_delta(value[1]);
        m_z.set_delta(value[2]);
    }

    Vector3D forces() const noexcept
    {
        return Vector3D(m_x.target(), m_y.target(), m_z.target());
    }

    void set_forces(Vector3D value) noexcept
    {
        m_x.set_target(value[0]);
        m_y.set_target(value[1]);
        m_z.set_target(value[2]);
    }

public:     // operators
    Parameter& operator[](const std::string& name) noexcept
    {
        if (name == "x") {
            return m_x;
        } else if (name == "y") {
            return m_y;
        } else if (name == "z") {
            return m_z;
        }

        return m_parameters[name];
    }

public:     // methods
    bool has_parameter(const std::string& name) noexcept
    {
        return m_parameters.find(name) != m_parameters.end();
    }
};

} // namespace EQlib