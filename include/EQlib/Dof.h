#pragma once

#include <vector>

namespace EQlib {

class Dof
{
private:    // variables
    double* m_ref_value;
    double* m_act_value;
    double* m_lower_bound;
    double* m_upper_bound;
    double* m_target;
    double* m_result;
    bool m_isfixed;

public:     // constructors
    Dof() noexcept
    : m_ref_value(nullptr)
    , m_act_value(nullptr)
    , m_lower_bound(nullptr)
    , m_upper_bound(nullptr)
    , m_target(nullptr)
    , m_result(nullptr)
    , m_isfixed(false)
    { }

    Dof(
        double* const ref_value,
        double* const act_value,
        double* const lower_bound,
        double* const upper_bound,
        double* const target,
        double* const result,
        const bool isfixed) noexcept
    : m_ref_value(ref_value)
    , m_act_value(act_value)
    , m_lower_bound(lower_bound)
    , m_upper_bound(upper_bound)
    , m_target(target)
    , m_result(result)
    , m_isfixed(isfixed)
    { }

public:     // getters and setters
    double ref_value() const
    {
        return *m_ref_value;
    }

    double act_value() const
    {
        return *m_act_value;
    }

    double lower_bound() const
    {
        return *m_lower_bound;
    }

    double upper_bound() const
    {
        return *m_upper_bound;
    }

    double delta() const
    {
        return *m_act_value - *m_ref_value;
    }

    void set_delta(double value) const
    {
        *m_act_value = *m_ref_value + value;
    }

    double residual() const
    {
        return *m_result - *m_target;
    }

    void set_residual(double value) const
    {
        *m_result = *m_target + value;
    }

    double target() const
    {
        return *m_target;
    }

    bool isfixed() const
    {
        return m_isfixed;
    }

    void set_isfixed(bool value)
    {
        m_isfixed = value;
    }

public:     // comparison
    bool operator==(const Dof& other) const noexcept
    {
        return m_act_value == other.m_act_value;
    }

    size_t hash() const noexcept
    {
        return (size_t)m_act_value;
    }

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Type = EQlib::Dof;

        py::class_<Type>(m, "Dof")
            // methods
            .def("__eq__", &Type::operator==)
            .def("__hash__", &Type::hash)
            // properties
            .def_property("delta", &Type::delta, &Type::set_delta)
            .def_property("residual", &Type::residual, &Type::set_residual)
            // read-only properties
            .def_property_readonly("isfixed", &Type::isfixed)
            .def_property_readonly("upper_bound", &Type::upper_bound)
            .def_property_readonly("lower_bound", &Type::lower_bound)
            .def_property_readonly("target", &Type::target)
        ;
    }
};

} // namespace EQlib

namespace std {

template <>
struct hash<EQlib::Dof>
{
    std::size_t operator()(const EQlib::Dof& dof) const noexcept
    {
        return dof.hash();
    }
};

} // namespace std