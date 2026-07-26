#pragma once

#include "Define.h"
#include "Variable.h"

#include <string>
#include <vector>

namespace eqlib {

class Objective {
private: // types
    using Type = Objective;

protected: // variables
    std::string m_name;
    bool m_is_active;
    std::vector<Pointer<Variable>> m_variables;

public: // constructors
    Objective()
        : m_is_active(true)
        , m_name("")
    {
    }

    Objective(const index nb_variables)
        : m_variables(nb_variables)
        , m_is_active(true)
        , m_name("")
    {
    }

    virtual ~Objective() = default;

public: // methods
    const Pointer<Variable>& variable(const index i) const
    {
        return m_variables[i];
    }

    const index nb_variables() const
    {
        return length(m_variables);
    }

    const std::vector<Pointer<Variable>>& variables() const
    {
        return m_variables;
    }

    virtual double compute(Ref<Vector> g, Ref<Matrix> h) const = 0;

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
    void set_variables(const std::vector<Pointer<Variable>>& value)
    {
        m_variables = value;
    }
};

} // namespace eqlib
