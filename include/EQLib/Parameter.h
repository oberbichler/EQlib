#pragma once

#include "Dof.h"

namespace EQLib {

class Parameter {
private:    // member variables
    double m_ref_value;
    double m_act_value;
    double m_target;
    double m_result;

public:     // constructors
    Parameter()
    : m_ref_value(0.0)
    , m_act_value(0.0)
    , m_target(0.0)
    , m_result(0.0)
    { }
    
    Parameter(
        const Parameter& other)
    : m_ref_value(other.m_ref_value)
    , m_act_value(other.m_act_value)
    , m_target(other.m_target)
    , m_result(other.m_result)
    { }
    
    Parameter(
        const Parameter&& other)
    : m_ref_value(std::move(other.m_ref_value))
    , m_act_value(std::move(other.m_act_value))
    , m_target(std::move(other.m_target))
    , m_result(std::move(other.m_result))
    { }

    Parameter(
        double ref_value,
        double act_value,
        double target,
        double result)
    : m_ref_value(ref_value)
    , m_act_value(act_value)
    , m_target(target)
    , m_result(result)
    { }

    Parameter(
        double value,
        double target)
    : Parameter(value, value, target, 0.0)
    { }

public:     // getters and setters
    double ref_value() const {
        return m_ref_value;
    }
    
    void set_ref_value(double value) {
        m_ref_value = value;
    }

    double act_value() const {
        return m_act_value;
    }
    
    void set_act_value(double value) {
        m_act_value = value;
    }

    double target() const {
        return m_target;
    }
    
    void set_target(double value) {
        m_target = value;
    }

    double result() const {
        return m_result;
    }
    
    void set_result(double value) {
        m_result = value;
    }

public:     // methods
    Dof dof() {
        return Dof(&m_ref_value, &m_act_value, &m_target, &m_result);
    }
};

} // namespace EQLib