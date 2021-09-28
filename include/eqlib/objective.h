#pragma once

#include "common.h"
#include "request.h"
#include "variable.h"

#include <string>
#include <vector>

namespace eqlib {

struct Objective {
    Objective()
        : m_is_active(true)
        , m_name("")
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

    // variables

    std::vector<Pointer<Variable>> m_variables;

    const index nb_variables() const
    {
        return len(m_variables);
    }

    const Pointer<Variable>& variable(const index i) const
    {
        return m_variables[i];
    }

    const std::vector<Pointer<Variable>>& variables() const
    {
        return m_variables;
    }

    void set_variables(std::vector<Pointer<Variable>>& value)
    {
        m_variables = value;
    }

    void add_variable(Pointer<Variable> value)
    {
        m_variables.push_back(value);
    }

    // variable_values

    Vector variable_values() const
    {
        Vector values(nb_variables());
        for (index i = 0; i < nb_variables(); i++) {
            values(i) = variable(i)->value();
        }
        return values;
    }

    void set_variable_values(const Ref<Vector> values)
    {
        assert(len(values) == nb_variables());

        for (index i = 0; i < nb_variables(); i++) {
            variable(i)->set_value(values(i));
        }
    }

    // compute

    virtual double compute(Ref<Vector> df, Ref<Matrix> hm, const Request request) = 0;
};

} // namespace eqlib
