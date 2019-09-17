#pragma once

#include "Define.h"

#include <limits>
#include <string>

namespace EQlib {

class Variable
{
private:    // static variables
    const static inline auto infinity = std::numeric_limits<double>::infinity();

private:    // variables
    double m_ref_value;
    double m_act_value;
    double m_lower_bound;
    double m_upper_bound;
    bool m_is_fixed;
    double m_multiplier;
    std::string m_name;

public:     // constructors
    Variable(
        const double ref_value,
        const double act_value,
        const double lower_bound,
        const double upper_bound,
        const bool is_fixed,
        const double multiplier,
        const std::string name) noexcept
    : m_ref_value(ref_value)
    , m_act_value(act_value)
    , m_lower_bound(lower_bound)
    , m_upper_bound(upper_bound)
    , m_is_fixed(is_fixed)
    , m_multiplier(multiplier)
    , m_name(name)
    { }

    Variable() noexcept
    : Variable(0.0, 0.0, -infinity, infinity, false, 1.0, "")
    { }

    Variable(
        const double value,
        const double lower_bound,
        const double upper_bound,
        const bool is_fixed,
        const std::string name) noexcept
    : Variable(value, value, lower_bound, upper_bound, is_fixed, 1.0, "")
    { }

    Variable(
        const double value) noexcept
    : Variable(value, value, -infinity, infinity, false, 1.0, "")
    { }

public:     // getters and setters
    double ref_value() const noexcept
    {
        return m_ref_value;
    }

    void set_ref_value(const double value) noexcept
    {
        m_ref_value = value;
    }

    double act_value() const noexcept
    {
        return m_act_value;
    }

    void set_act_value(const double value) noexcept
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

    double delta() const noexcept
    {
        return m_act_value - m_ref_value;
    }

    void set_delta(const double value) noexcept
    {
        m_act_value = m_ref_value + value;
    }

    bool is_fixed() const noexcept
    {
        return m_is_fixed;
    }

    void set_is_fixed(const bool value) noexcept
    {
        m_is_fixed = value;
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

public:     // methods
    std::string to_string() const noexcept
    {
        if (m_name.empty()) {
            return format(
                "<Variable value={} is_fixed={} bounds=({}, {}) at {:#x}>",
                act_value(), is_fixed(), lower_bound(), upper_bound(),
                size_t(this));
        } else {
            return format(
                "<Variable '{}' value={} is_fixed={} bounds=({}, {}) at {:#x}>",
                name(), act_value(), is_fixed(), lower_bound(), upper_bound(),
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
        return act_value();
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
            .def(py::init<double, double, double, double, bool, double, std::string>(),
                "ref_value"_a, "act_value"_a, "lower_bound"_a=-infinity,
                "upper_bound"_a=infinity, "is_fixed"_a=false,
                "multiplier"_a=0.0, "name"_a="")
            .def(py::init<double, double, double, bool, std::string>(),
                "value"_a, "lower_bound"_a=-infinity, "upper_bound"_a=infinity,
                "is_fixed"_a=false, "name"_a="")
            .def(py::init<>())
            .def_property("ref_value", &Type::ref_value, &Type::set_ref_value)
            .def_property("act_value", &Type::act_value, &Type::set_act_value)
            .def_property("lower_bound", &Type::lower_bound,
                &Type::set_lower_bound)
            .def_property("upper_bound", &Type::upper_bound,
                &Type::set_upper_bound)
            .def_property("delta", &Type::delta, &Type::set_delta)
            .def_property("is_fixed", &Type::is_fixed, &Type::set_is_fixed)
            .def_property("multiplier", &Type::multiplier,
                &Type::set_multiplier)
            .def_property("name", &Type::name, &Type::set_name)
        ;
    }
};

} // namespace EQlib