#pragma once

#include <vector>

namespace EQLib {

class Dof {
private:    // member variables
    double* m_ref_value;
    double* m_act_value;
    double* m_target;
    double* m_result;

public:     // constructors
    Dof()
    : m_ref_value(nullptr)
    , m_act_value(nullptr)
    , m_target(nullptr)
    , m_result(nullptr)
    { }
    
    Dof(
        double* const ref_value,
        double* const act_value,
        double* const target,
        double* const result)
    : m_ref_value(ref_value)
    , m_act_value(act_value)
    , m_target(target)
    , m_result(result)
    { }

public:     // getters and setters
    double delta() const {
        return *m_act_value - *m_ref_value;
    }
    
    void set_delta(double value) const {
        *m_act_value = *m_ref_value + value;
    }
    
    double residual() const {
        return *m_result - *m_target;
    }
    
    void set_residual(double value) const {
        *m_result = *m_target + value;
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

} // namespace EQLib