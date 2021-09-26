#pragma once

#include "common.h"
#include "equation.h"
#include "variable.h"

#include <vector>
#include <string>

namespace eqlib
{

    struct Constraint
    {
        Constraint() : m_is_active(true), m_name("")
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

        const Pointer<Equation> &equation(const index i) const
        {
            return m_equations[i];
        }

        const std::vector<Pointer<Equation>> &equations() const
        {
            return m_equations;
        }

        const void set_equations(std::vector<Pointer<Equation>> &value)
        {
            m_equations = value;
        }

        // variables

        std::vector<Pointer<Variable>> m_variables;

        const index nb_variables() const
        {
            return len(m_variables);
        }

        const Pointer<Variable> &variable(const index i) const
        {
            return m_variables[i];
        }

        const std::vector<Pointer<Variable>> &variables() const
        {
            return m_variables;
        }

        const void set_variables(std::vector<Pointer<Variable>> &value)
        {
            m_variables = value;
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

        virtual void compute(const Ref<Vector> g, const Ref<Matrix> dg, const Ref<Matrix> hm, const Request request) = 0;
    };

} // namespace eqlib
