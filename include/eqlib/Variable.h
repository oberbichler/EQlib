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
        const std::string name) noexcept;

    Variable() noexcept;

    Variable(
        const double value,
        const double lower_bound,
        const double upper_bound,
        const bool is_active,
        const std::string name) noexcept;

    Variable(
        const double value) noexcept;

public:     // methods
    double value() const noexcept;

    void set_value(const double value) noexcept;

    double lower_bound() const noexcept;

    void set_lower_bound(const double value) noexcept;

    double upper_bound() const noexcept;

    void set_upper_bound(const double value) noexcept;

    bool is_active() const noexcept;

    void set_active(const bool value) noexcept;

    double multiplier() const noexcept;

    void set_multiplier(const double value) noexcept;

    const std::string& name() const noexcept

    void set_name(const std::string& value) noexcept;

    std::string to_string() const noexcept;

public:     // comparison
    bool operator==(const Variable& other) const noexcept;

    size_t hash() const noexcept;

public:     // operators
    operator double();

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