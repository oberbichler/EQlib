#pragma once

#include "Define.h"
#include "Dof.h"

#include <limits>
#include <string>

namespace EQlib {

class Parameter
{
private:    // variables
    double m_ref_value;
    double m_act_value;
    double m_min_value;
    double m_max_value;
    double m_target;
    double m_result;
    bool m_isfixed;
    std::string m_name;

public:     // constructors
    Parameter(
        double ref_value,
        double act_value,
        double target,
        double result,
        bool isfixed) noexcept
    : m_ref_value(ref_value)
    , m_act_value(act_value)
    , m_min_value(-std::numeric_limits<double>::infinity())
    , m_max_value(std::numeric_limits<double>::infinity())
    , m_target(target)
    , m_result(result)
    , m_isfixed(isfixed)
    { }

    Parameter() noexcept
    : Parameter(0, 0, 0, 0, false)
    { }

    Parameter(
        double value,
        double target=0,
        bool isfixed=false) noexcept
    : Parameter(value, value, target, 0, isfixed)
    { }

public:     // getters and setters
    double ref_value() const noexcept
    {
        return m_ref_value;
    }

    void set_ref_value(double value) noexcept
    {
        m_ref_value = value;
    }

    double act_value() const noexcept
    {
        return m_act_value;
    }

    void set_act_value(double value) noexcept
    {
        m_act_value = value;
    }

    double min_value() const noexcept
    {
        return m_min_value;
    }

    void set_min_value(double value) noexcept
    {
        m_min_value = value;
    }

    double max_value() const noexcept
    {
        return m_max_value;
    }

    void set_max_value(double value) noexcept
    {
        m_max_value = value;
    }

    double delta() const noexcept
    {
        return m_act_value - m_ref_value;
    }

    void set_delta(double value) noexcept
    {
        m_act_value = m_ref_value + value;
    }

    double target() const noexcept
    {
        return m_target;
    }

    void set_target(double value) noexcept
    {
        m_target = value;
    }

    double result() const noexcept
    {
        return m_result;
    }

    void set_result(double value) noexcept
    {
        m_result = value;
    }

    double residual() const noexcept
    {
        return m_target - m_result;
    }

    void set_residual(double value) noexcept
    {
        m_result = m_target - value;
    }

    bool isfixed() const noexcept
    {
        return m_isfixed;
    }

    void set_isfixed(bool value) noexcept
    {
        m_isfixed = value;
    }

    std::string name() const noexcept
    {
        return m_name;
    }

    void set_name(const std::string& value) noexcept
    {
        m_name = value;
    }

public:     // methods
    Dof dof() noexcept
    {
        return Dof(&m_ref_value, &m_act_value, &m_min_value, &m_max_value,
            &m_target, &m_result, m_isfixed);
    }

    std::string to_string() const noexcept
    {
        if (m_name.empty()) {
            return format(
                "<Parameter value={} isfixed={} bounds=({}, {}) at {:#x}>",
                act_value(), isfixed(), min_value(), max_value(),
                size_t(&m_act_value));
        } else {
            return format(
                "<Parameter '{}' value={} isfixed={} bounds=({}, {}) at {:#x}>",
                name(), act_value(), isfixed(), min_value(), max_value(),
                size_t(&m_act_value));
        }
    }
};

} // namespace EQlib