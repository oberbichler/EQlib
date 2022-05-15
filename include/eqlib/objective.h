#pragma once

#include "common.h"
#include "request.h"
#include "parameter.h"

#include <string>
#include <vector>

namespace eqlib {

struct Objective {
    Objective()
        : m_is_active(true)
        , m_name("")
    {
    }

    Objective(const index nb_parameters)
        : m_is_active(true)
        , m_name("")
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

    // parameters

    Parameters m_parameters;

    const index nb_parameters() const
    {
        return len(m_parameters);
    }

    const Pointer<Parameter>& parameter(const index i) const
    {
        return m_parameters[i];
    }

    const Parameters& parameters() const
    {
        return m_parameters;
    }

    void set_parameters(Parameters& value)
    {
        m_parameters = value;
    }

    void add_parameter(Pointer<Parameter> value)
    {
        m_parameters.push_back(value);
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

    virtual double compute(Ref<Vector> df, Ref<Matrix> hm, const Request request) = 0;
};

} // namespace eqlib
