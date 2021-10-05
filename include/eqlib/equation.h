#pragma once

#include "common.h"

#include <string>
#include <vector>

namespace eqlib {

struct Equation {
    Equation()
        : m_value(0)
        , m_is_active(true)
        , m_name("")
        , m_lower_bound(LOWEST)
        , m_upper_bound(HIGHEST)
    {
    }

    Equation(const bool is_active, const std::string& name)
        : m_value(0)
        , m_is_active(is_active)
        , m_name(name)
        , m_lower_bound(LOWEST)
        , m_upper_bound(HIGHEST)
    {
    }

    // value

    double m_value;

    double value() const
    {
        return m_value;
    }

    void set_value(const double value)
    {
        m_value = value;
    }

    // is_active

    bool m_is_active;

    bool is_active() const
    {
        return m_is_active;
    }

    void set_is_active(const bool is_active)
    {
        m_is_active = is_active;
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

    // lower bound

    double m_lower_bound;

    double lower_bound() const
    {
        return m_lower_bound;
    }

    void set_lower_bound(const double value)
    {
        m_lower_bound = value;
    }

    // upper bound

    double m_upper_bound;

    double upper_bound() const
    {
        return m_upper_bound;
    }

    void set_upper_bound(const double value)
    {
        m_upper_bound = value;
    }
};

using Equations = std::vector<Pointer<Equation>>;

} // namespace eqlib
