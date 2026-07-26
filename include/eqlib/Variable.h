#pragma once

#include "Define.h"

#include <string>

namespace eqlib {

class Variable {
private: // types
    using Type = Variable;

private: // variables
    double m_act_value;
    double m_lower_bound;
    double m_upper_bound;
    bool m_is_active;
    double m_multiplier;
    std::string m_name;

public: // constructors
    Variable(
        const double value,
        const double lower_bound,
        const double upper_bound,
        const bool is_active,
        const double multiplier,
        const std::string name) noexcept
        : m_act_value(value)
        , m_lower_bound(lower_bound)
        , m_upper_bound(upper_bound)
        , m_is_active(is_active)
        , m_multiplier(multiplier)
        , m_name(name)
    {
    }

    Variable() noexcept
        : Variable(0.0, -infinity, infinity, true, 1.0, "")
    {
    }

    Variable(
        const double value,
        const double lower_bound,
        const double upper_bound,
        const bool is_active,
        const std::string name) noexcept
        : Variable(value, lower_bound, upper_bound, is_active, 1.0, name)
    {
    }

    Variable(
        const double value) noexcept
        : Variable(value, -infinity, infinity, true, 1.0, "")
    {
    }

public: // methods
    double value() const noexcept
    {
        return m_act_value;
    }

    double& value() noexcept
    {
        return m_act_value;
    }

    void set_value(const double value) noexcept
    {
        m_act_value = value;
    }

    double lower_bound() const noexcept
    {
        return m_lower_bound;
    }

    void set_lower_bound(const double value) noexcept
    {
        m_lower_bound = value;
    }

    double upper_bound() const noexcept
    {
        return m_upper_bound;
    }

    void set_upper_bound(const double value) noexcept
    {
        m_upper_bound = value;
    }

    bool is_active() const noexcept
    {
        return m_is_active;
    }

    void set_active(const bool value) noexcept
    {
        m_is_active = value;
    }

    double multiplier() const noexcept
    {
        return m_multiplier;
    }

    void set_multiplier(const double value) noexcept
    {
        m_multiplier = value;
    }

    void clamp() noexcept
    {
        if (value() < lower_bound()) {
            set_value(lower_bound());
        } else if (upper_bound() < value()) {
            set_value(upper_bound());
        }
    }

    const std::string& name() const noexcept
    {
        return m_name;
    }

    void set_name(const std::string& value) noexcept
    {
        m_name = value;
    }

    std::string to_string() const noexcept
    {
        if (m_name.empty()) {
            return format(
                "<Variable value={} is_active={} bounds=({}, {}) at {:#x}>",
                value(), is_active(), lower_bound(), upper_bound(),
                size_t(this));
        } else {
            return format(
                "<Variable '{}' value={} is_active={} bounds=({}, {}) at {:#x}>",
                name(), value(), is_active(), lower_bound(), upper_bound(),
                size_t(this));
        }
    }

public: // comparison
    bool operator==(const Variable& other) const noexcept
    {
        return this == &other;
    }

    size_t hash() const noexcept
    {
        return (size_t)this;
    }

public: // operators
    operator double()
    {
        return value();
    }
};

} // namespace eqlib
