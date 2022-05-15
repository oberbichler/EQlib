#pragma once

#include "parameter.h"

namespace eqlib {

struct Node {
    Node()
        : m_name("")
    {
        m_x = new_<Parameter>(0);
        m_y = new_<Parameter>(0);
        m_z = new_<Parameter>(0);
    }

    Node(const double x, const double y, const double z, const std::string& name)
        : m_name(name)
    {
        m_x = new_<Parameter>(x);
        m_y = new_<Parameter>(y);
        m_z = new_<Parameter>(z);
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

    Pointer<Parameter> m_x;

    Pointer<Parameter> x()
    {
        return m_x;
    }

    // y

    Pointer<Parameter> m_y;

    Pointer<Parameter> y()
    {
        return m_y;
    }

    // z

    Pointer<Parameter> m_z;

    Pointer<Parameter> z()
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
