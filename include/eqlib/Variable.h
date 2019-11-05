#pragma once

#include "Define.h"

#include <limits>
#include <string>

namespace eqlib {

class Variable
{
private:    // static variables
    const static inline auto infinity = std::numeric_limits<double>::infinity();

private:    // variables
    double m_act_value;
    double m_lower_bound;
    double m_upper_bound;
    bool m_is_active;
    double m_multiplier;
    std::string m_name;

public:     // constructors
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
    { }

    Variable() noexcept
    : Variable(0.0, -infinity, infinity, true, 1.0, "")
    { }

    Variable(
        const double value,
        const double lower_bound,
        const double upper_bound,
        const bool is_active,
        const std::string name) noexcept
    : Variable(value, lower_bound, upper_bound, is_active, 1.0, name)
    { }

    Variable(
        const double value) noexcept
    : Variable(value, -infinity, infinity, true, 1.0, "")
    { }

public:     // methods
    double value() const noexcept
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

    std::string name() const noexcept
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

public:     // comparison
    bool operator==(const Variable& other) const noexcept
    {
        return this == &other;
    }

    size_t hash() const noexcept
    {
        return (size_t)this;
    }

public:     // operators
    operator double()
    {
        return value();
    }

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Type = Variable;
        using Holder = Pointer<Type>;

        py::class_<Type, Holder>(m, "Variable")
            .def(py::init<double, double, double, bool, double, std::string>(),
                "value"_a, "lower_bound"_a=-infinity,
                "upper_bound"_a=infinity, "is_active"_a=true,
                "multiplier"_a=1.0, "name"_a="")
            .def(py::init<>())
            .def_property("value", &Type::value, &Type::set_value)
            .def_property("lower_bound", &Type::lower_bound,
                &Type::set_lower_bound)
            .def_property("upper_bound", &Type::upper_bound,
                &Type::set_upper_bound)
            .def_property("is_active", &Type::is_active, &Type::set_active)
            .def_property("multiplier", &Type::multiplier,
                &Type::set_multiplier)
            .def_property("name", &Type::name, &Type::set_name)
            .def("__float__", &Type::operator double)
        ;
    }
};

} // namespace eqlib