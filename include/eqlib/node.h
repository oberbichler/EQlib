#pragma once

#include "variable.h"

namespace eqlib {

struct Node {
    Node()
        : m_name("")
    {
        m_x = new_<Variable>(0);
        m_y = new_<Variable>(0);
        m_z = new_<Variable>(0);
    }

    Node(const double x, const double y, const double z, const std::string& name)
        : m_name(name)
    {
        m_x = new_<Variable>(x);
        m_y = new_<Variable>(y);
        m_z = new_<Variable>(z);
    }

    // name

    std::string m_name;

    std::string name() const
    {
        return m_name;
    }

    void set_name(const std::string value)
    {
        m_name = value;
    }

    // x

    Pointer<Variable> m_x;

    Pointer<Variable> x()
    {
        return m_x;
    }

    // y

    Pointer<Variable> m_y;

    Pointer<Variable> y()
    {
        return m_y;
    }

    // z

    Pointer<Variable> m_z;

    Pointer<Variable> z()
    {
        return m_z;
    }

    // location

    Vector3 location()
    {
        return Vector3(m_x->value(), m_y->value(), m_z->value());
    }

    void set_location(const Vector3 value)
    {
        m_x->set_value(value(0));
        m_y->set_value(value(1));
        m_z->set_value(value(2));
    }
};

} // namespace eqlib
