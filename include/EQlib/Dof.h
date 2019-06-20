#pragma once

#include <vector>

namespace EQlib {

class Dof
{
private:    // variables
    double* m_ref_value;
    double* m_act_value;
    double* m_min_value;
    double* m_max_value;
    double* m_target;
    double* m_result;
    bool m_isfixed;

public:     // constructors
    Dof() noexcept
    : m_ref_value(nullptr)
    , m_act_value(nullptr)
    , m_min_value(nullptr)
    , m_max_value(nullptr)
    , m_target(nullptr)
    , m_result(nullptr)
    , m_isfixed(false)
    { }

    Dof(
        double* const ref_value,
        double* const act_value,
        double* const min_value,
        double* const max_value,
        double* const target,
        double* const result,
        const bool isfixed) noexcept
    : m_ref_value(ref_value)
    , m_act_value(act_value)
    , m_min_value(min_value)
    , m_max_value(max_value)
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

    double min_value() const
    {
        return *m_min_value;
    }

    double max_value() const
    {
        return *m_max_value;
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