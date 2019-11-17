#pragma once

#include "Define.h"

#include <limits>
#include <string>

namespace eqlib {

Variable::Variable(
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

Variable::Variable() noexcept : Variable(0.0, -infinity, infinity, true, 1.0, "")
{   
}

Variable::Variable(
    const double value,
    const double lower_bound,
    const double upper_bound,
    const bool is_active,
    const std::string name) noexcept
: Variable(value, lower_bound, upper_bound, is_active, 1.0, name)
{
}

Variable::Variable(const double value) noexcept
: Variable(value, -infinity, infinity, true, 1.0, "")
{
}

double Variable::value() const noexcept
{
    return m_act_value;
}

void Variable::set_value(const double value) noexcept
{
    m_act_value = value;
}

double Variable::lower_bound() const noexcept
{
    return m_lower_bound;
}

void Variable::set_lower_bound(const double value) noexcept
{
    m_lower_bound = value;
}

double Variable::upper_bound() const noexcept
{
    return m_upper_bound;
}

void Variable::set_upper_bound(const double value) noexcept
{
    m_upper_bound = value;
}

bool Variable::is_active() const noexcept
{
    return m_is_active;
}

void Variable::set_active(const bool value) noexcept
{
    m_is_active = value;
}

double Variable::multiplier() const noexcept
{
    return m_multiplier;
}

void Variable::set_multiplier(const double value) noexcept
{
    m_multiplier = value;
}

std::string name() const noexcept
{
    return m_name;
}

void Variable::set_name(const std::string& value) noexcept
{
    m_name = value;
}

std::string Variable::to_string() const noexcept
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

bool Variable::operator==(const Variable& other) const noexcept
{
    return this == &other;
}

size_t Variable::hash() const noexcept
{
    return (size_t)this;
}

Variable::operator double()
{
    return value();
}

} // namespace eqlib