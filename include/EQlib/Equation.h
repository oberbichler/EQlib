#pragma once

#include "Define.h"

#include <limits>
#include <string>

namespace EQlib {

class Equation
{
private:    // static variables
    const static inline auto infinity = std::numeric_limits<double>::infinity();

private:    // variables
    double m_lower_bound;
    double m_upper_bound;
    bool m_is_active;
    double m_multiplier;
    std::string m_name;

public:     // constructors
    Equation(
        const double lower_bound,
        const double upper_bound,
        const double multiplier,
        const std::string name) noexcept
    : m_lower_bound(lower_bound)
    , m_upper_bound(upper_bound)
    , m_is_active(true)
    , m_multiplier(multiplier)
    , m_name(name)
    { }

    Equation() noexcept
    : Equation(-infinity, infinity, 0.0, "")
    { }

public:     // getters and setters
    bool is_active() const noexcept
    {
        return m_is_active;
    }

    void set_is_active(const bool value) noexcept
    {
        m_is_active = value;
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
                "<Equation bounds=({}, {}) at {:#x}>", lower_bound(),
                    upper_bound(), size_t(this));
        } else {
            return format(
                "<Equation '{}' bounds=({}, {}) at {:#x}>", name(),
                    lower_bound(), upper_bound(), size_t(this));
        }
    }

public:     // comparison
    bool operator==(const Equation& other) const noexcept
    {
        return this == &other;
    }

    size_t hash() const noexcept
    {
        return (size_t)this;
    }

public:     // operators

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Type = Equation;
        using Holder = Pointer<Type>;

        py::class_<Type, Holder>(m, "Equation")
            .def(py::init<double, double, double, std::string>(),
                "lower_bound"_a=-infinity, "upper_bound"_a=infinity,
                "multiplier"_a=0.0, "name"_a="")
            .def(py::init<>())
            .def_property("is_active", &Type::is_active, &Type::set_is_active)
            .def_property("lower_bound", &Type::lower_bound,
                &Type::set_lower_bound)
            .def_property("upper_bound", &Type::upper_bound,
                &Type::set_upper_bound)
            .def_property("name", &Type::name, &Type::set_name)
            .def_property("multiplier", &Type::multiplier,
                &Type::set_multiplier)
            .def("__repr__", &Type::to_string)
        ;
    }
};

} // namespace EQlib