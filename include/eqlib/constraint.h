#pragma once

#include "common.h"
#include "equation.h"
#include "request.h"
#include "parameter.h"

#include <string>
#include <vector>

namespace eqlib {

struct Constraint {
    Constraint()
        : m_is_active(true)
        , m_name("")
    {
    }

    Constraint(const index nb_equations, const index nb_parameters)
        : m_is_active(true)
        , m_name("")
        , m_equations(nb_equations)
        , m_parameters(nb_parameters)
    {
    }

    // is_active

    bool m_is_active;

    bool is_active() const
    {
        return m_is_active;
    }

    void set_is_active(const bool is_active)
    {
        m_is_active = is_active;
    }

    // name

    std::string m_name;

    std::string name() const
    {
        return m_name;
    }

    void set_name(const std::string value)
    {
        m_name = value;
    }

    // equations

    std::vector<Pointer<Equation>> m_equations;

    const index nb_equations() const
    {
        return len(m_equations);
    }

    const Pointer<Equation>& equation(const index i) const
    {
        return m_equations[i];
    }

    const std::vector<Pointer<Equation>>& equations() const
    {
        return m_equations;
    }

    const void set_equations(std::vector<Pointer<Equation>>& value)
    {
        m_equations = value;
    }

    // parameters

    std::vector<Pointer<Parameter>> m_parameters;

    const index nb_parameters() const
    {
        return len(m_parameters);
    }

    const Pointer<Parameter>& parameter(const index i) const
    {
        return m_parameters[i];
    }

    const std::vector<Pointer<Parameter>>& parameters() const
    {
        return m_parameters;
    }

    const void set_parameters(std::vector<Pointer<Parameter>>& value)
    {
        m_parameters = value;
    }

    // parameter_values

    Vector parameter_values() const
    {
        Vector values(nb_parameters());
        for (index i = 0; i < nb_parameters(); i++) {
            values(i) = parameter(i)->value();
        }
        return values;
    }

    void set_parameter_values(const Ref<Vector> values)
    {
        assert(len(values) == nb_parameters());

        for (index i = 0; i < nb_parameters(); i++) {
            parameter(i)->set_value(values(i));
        }
    }

    // compute

    virtual void compute(const Ref<Vector> g, const Ref<Matrix> dg, const Ref<Matrix> hm, const Request request) = 0;
};

} // namespace eqlib
