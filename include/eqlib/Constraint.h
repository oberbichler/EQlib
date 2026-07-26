#pragma once

#include "Define.h"
#include "Equation.h"
#include "Variable.h"

#include <string>
#include <vector>

namespace eqlib {

class Constraint {
private: // types
    using Type = Constraint;

protected: // variables
    std::string m_name;
    bool m_is_active;
    std::vector<Pointer<Equation>> m_equations;
    std::vector<Pointer<Variable>> m_variables;

public: // constructors
    Constraint()
        : m_is_active(true)
        , m_name("")
    {
    }

    Constraint(const index nb_equations, const index nb_variables)
        : m_equations(nb_equations)
        , m_variables(nb_variables)
        , m_is_active(true)
        , m_name("")
    {
    }

    virtual ~Constraint() = default;

public: // methods
    const Pointer<Equation>& equation(const index i) const
    {
        return m_equations[i];
    }

    const Pointer<Variable>& variable(const index i) const
    {
        return m_variables[i];
    }

    const index nb_equations() const
    {
        return length(m_equations);
    }

    const index nb_variables() const
    {
        return length(m_variables);
    }

    const std::vector<Pointer<Equation>>& equations() const
    {
        return m_equations;
    }

    const std::vector<Pointer<Variable>>& variables() const
    {
        return m_variables;
    }

    virtual void compute(Ref<Vector> rs, const std::vector<Ref<Vector>>& gs, const std::vector<Ref<Matrix>>& hs) const = 0;

    bool is_active() const noexcept
    {
        return m_is_active;
    }

    void set_active(const bool value) noexcept
    {
        m_is_active = value;
    }

    const std::string& name() const
    {
        return m_name;
    }

    void set_name(const std::string& value)
    {
        m_name = value;
    }

public: // methods: setters
    void set_equations(const std::vector<Pointer<Equation>>& value)
    {
        m_equations = value;
    }

    void set_variables(const std::vector<Pointer<Variable>>& value)
    {
        m_variables = value;
    }
};

} // namespace eqlib
