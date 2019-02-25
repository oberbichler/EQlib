#pragma once

#include "Dof.h"

namespace EQlib {

class Parameter {
private:    // member variables
    double m_ref_value;
    double m_act_value;
    double m_target;
    double m_result;
    bool m_isfixed;

public:     // constructors
    Parameter(
        double ref_value,
        double act_value,
        double target,
        double result,
        bool isfixed)
    : m_ref_value(ref_value)
    , m_act_value(act_value)
    , m_target(target)
    , m_result(result)
    , m_isfixed(isfixed)
    { }

    Parameter()
    : Parameter(0, 0, 0, 0, false)
    { }

    Parameter(
        double value,
        double target=0,
        bool isfixed=false)
    : Parameter(value, value, target, 0, isfixed)
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

    double delta() const {
        return m_act_value - m_ref_value;
    }

    void set_delta(double value) {
        m_act_value = m_ref_value + value;
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

    double residual() const {
        return m_target - m_result;
    }

    void set_residual(double value) {
        m_result = m_target - value;
    }

    bool isfixed() const {
        return m_isfixed;
    }

    void set_isfixed(bool value) {
        m_isfixed = value;
    }

public:     // methods
    Dof dof() {
        return Dof(&m_ref_value, &m_act_value, &m_target, &m_result, m_isfixed);
    }
};

} // namespace EQlib