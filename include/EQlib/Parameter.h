#pragma once

#include "Define.h"
#include "Parameter.h"

#include <limits>
#include <string>

namespace EQlib {

class Parameter
{
private:    // variables
    double m_ref_value;
    double m_act_value;
    double m_lower_bound;
    double m_upper_bound;
    double m_target;
    double m_result;
    bool m_is_fixed;
    std::string m_name;

public:     // constructors
    Parameter(
        const double ref_value,
        const double act_value,
        const double target,
        const double result,
        const bool is_fixed,
        const std::string name) noexcept
    : m_ref_value(ref_value)
    , m_act_value(act_value)
    , m_lower_bound(-std::numeric_limits<double>::infinity())
    , m_upper_bound(std::numeric_limits<double>::infinity())
    , m_target(target)
    , m_result(result)
    , m_is_fixed(is_fixed)
    , m_name(name)
    { }

    Parameter() noexcept
    : Parameter(0, 0, 0, 0, false, "")
    { }

    Parameter(
        const double value,
        const double target=0,
        const bool is_fixed=false) noexcept
    : Parameter(value, value, target, 0, is_fixed, "")
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

    double target() const noexcept
    {
        return m_target;
    }

    void set_target(const double value) noexcept
    {
        m_target = value;
    }

    double result() const noexcept
    {
        return m_result;
    }

    void set_result(const double value) noexcept
    {
        m_result = value;
    }

    double residual() const noexcept
    {
        return m_target - m_result;
    }

    void set_residual(const double value) noexcept
    {
        m_result = m_target - value;
    }

    bool is_fixed() const noexcept
    {
        return m_is_fixed;
    }

    void set_is_fixed(const bool value) noexcept
    {
        m_is_fixed = value;
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
                "<Parameter value={} is_fixed={} bounds=({}, {}) at {:#x}>",
                act_value(), is_fixed(), lower_bound(), upper_bound(),
                size_t(this));
        } else {
            return format(
                "<Parameter '{}' value={} is_fixed={} bounds=({}, {}) at {:#x}>",
                name(), act_value(), is_fixed(), lower_bound(), upper_bound(),
                size_t(this));
        }
    }

public:     // comparison
    bool operator==(const Parameter& other) const noexcept
    {
        return &m_act_value == &(other.m_act_value);
    }

    size_t hash() const noexcept
    {
        return (size_t)&m_act_value;
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

        using Type = Parameter;
        using Holder = Pointer<Type>;

        py::class_<Type, Holder>(m, "Parameter")
            .def(py::init<double, double, double, double, bool, std::string>(),
                "ref_value"_a, "act_value"_a, "target"_a=0, "result"_a=0,
                "is_fixed"_a=false, "name"_a="")
            .def(py::init<double, double, bool>(), "value"_a, "target"_a=0,
                "is_fixed"_a=false)
            .def(py::init<>())
            .def_property("ref_value", &Type::ref_value, &Type::set_ref_value)
            .def_property("act_value", &Type::act_value, &Type::set_act_value)
            .def_property("lower_bound", &Type::lower_bound,
                &Type::set_lower_bound)
            .def_property("upper_bound", &Type::upper_bound,
                &Type::set_upper_bound)
            .def_property("delta", &Type::delta, &Type::set_delta)
            .def_property("target", &Type::target, &Type::set_target)
            .def_property("result", &Type::result, &Type::set_result)
            .def_property("residual", &Type::residual, &Type::set_residual)
            .def_property("is_fixed", &Type::is_fixed, &Type::set_is_fixed)
            .def_property("name", &Type::name, &Type::set_name)
            .def(py::pickle([](const Type& self) {
                    return py::make_tuple(
                        self.ref_value(),
                        self.act_value(),
                        self.target(),
                        self.result(),
                        self.is_fixed(),
                        self.lower_bound(),
                        self.upper_bound(),
                        self.name());
                }, [](py::tuple tuple) {
                    if (tuple.size() != 8) {
                        throw std::runtime_error("Invalid state!");
                    }

                    const auto ref_value = tuple[0].cast<double>();
                    const auto act_value = tuple[1].cast<double>();
                    const auto target = tuple[2].cast<double>();
                    const auto result = tuple[3].cast<double>();
                    const auto is_fixed = tuple[4].cast<bool>();
                    const auto lower_bound = tuple[5].cast<double>();
                    const auto upper_bound = tuple[6].cast<double>();
                    const auto name = tuple[7].cast<std::string>();

                    return Type(ref_value, act_value, target, result, is_fixed);
                }
            ))
            .def("__float__", [](const Type& self) { return self.act_value(); })
            .def("__repr__", &Type::to_string)
            .def("__copy__", [](const Type& self) { return Type(self); })
            .def("__deepcopy__", [](const Type& self, py::dict& memo) {
                return Type(self); }, "memodict"_a)
        ;
    }
};

} // namespace EQlib