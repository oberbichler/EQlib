#pragma once

#include <eqlib/Equation.h>

#include <limits>
#include <string>

namespace eqlib {

Equation::Equation(const double lower_bound, const double upper_bound, const double multiplier, const std::string name) noexcept
: m_lower_bound(lower_bound)
, m_upper_bound(upper_bound)
, m_is_active(true)
, m_multiplier(multiplier)
, m_name(name)
{ }

Equation::Equation() noexcept : Equation(-infinity, infinity, 0.0, "")
{
}

bool Equation::is_active() const noexcept
{
    return m_is_active;
}

void Equation::set_active(const bool value) noexcept
{
    m_is_active = value;
}

double Equation::lower_bound() const noexcept
{
    return m_lower_bound;
}

void Equation::set_lower_bound(const double value) noexcept
{
    m_lower_bound = value;
}

double Equation::upper_bound() const noexcept
{
    return m_upper_bound;
}

void Equation::set_upper_bound(const double value) noexcept
{
    m_upper_bound = value;
}

double Equation::multiplier() const noexcept
{
    return m_multiplier;
}

void Equation::set_multiplier(const double value) noexcept
{
    m_multiplier = value;
}

const std::string& Equation::name() const noexcept
{
    return m_name;
}

void Equation::set_name(const std::string& value) noexcept
{
    m_name = value;
}

std::string Equation::to_string() const noexcept
{
    if (m_name.empty()) {
        return format(
            "<Equation bounds=({}, {}) at {:#x}>", lower_bound(),
                upper_bound(), size_t(this));
    } else {
        return format(
            "<Equation '{}' bounds=({}, {}) at {:#x}>", name(),
                lower_bound(), upper_bound(), size_t(this));
    }
}

bool Equation::operator==(const Equation& other) const noexcept
{
    return this == &other;
}

size_t Equation::hash() const noexcept
{
    return (size_t)this;
}

} // namespace eqlib