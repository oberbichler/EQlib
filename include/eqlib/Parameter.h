#pragma once

#include "Define.h"

#include <string>

namespace eqlib {

class Parameter {
private: // types
    using Type = Parameter;

private: // variables
    double m_value;
    std::string m_name;

public: // constructors
    Parameter(
        const double value,
        const std::string name) noexcept
        : m_value(value)
        , m_name(name)
    {
    }

    Parameter() noexcept
        : Parameter(0.0, "")
    {
    }

    Parameter(const double value) noexcept
        : Parameter(value, "")
    {
    }

public: // methods
    double value() const noexcept
    {
        return m_value;
    }

    void set_value(const double value) noexcept
    {
        m_value = value;
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
            return format("<Parameter value={} at {:#x}>", value(), size_t(this));
        } else {
            return format("<Parameter '{}' value={} at {:#x}>", name(), value(), size_t(this));
        }
    }

public: // comparison
    bool operator==(const Parameter& other) const noexcept
    {
        return this == &other;
    }

    size_t hash() const noexcept
    {
        return (size_t)this;
    }

public: // operators
    operator double() const
    {
        return value();
    }

public: // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Holder = Pointer<Type>;

        py::class_<Type, Holder>(m, "Parameter")
            .def(py::init<>())
            .def(py::init<double, std::string>(), "value"_a = 0.0, "name"_a = "")
            .def_property("value", &Type::value, &Type::set_value)
            .def_property("name", &Type::name, &Type::set_name)
            .def("__float__", &Type::operator double);
    }
};

} // namespace eqlib