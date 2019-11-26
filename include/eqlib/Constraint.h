#pragma once

#include "Define.h"
#include "Equation.h"
#include "Variable.h"

#include <string>
#include <vector>

namespace eqlib {

class Constraint
{
protected:  // variables
    std::string m_name;
    bool m_is_active;
    std::vector<Pointer<Equation>> m_equations;
    std::vector<Pointer<Variable>> m_variables;

public:     // constructors
    Constraint() : m_is_active(true), m_name("")
    {
    }

    Constraint(const index nb_equations, const index nb_variables)
    : m_equations(nb_equations), m_variables(nb_variables), m_is_active(true), m_name("")
    {
    }

    virtual ~Constraint() = default;

public:     // methods
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

protected:  // methods
    void set_equations(const std::vector<Pointer<Equation>>& value)
    {
        m_equations = value;
    }

    void set_variables(const std::vector<Pointer<Variable>>& value)
    {
        m_variables = value;
    }

public:     // python
    template <typename T>
    class PyConstraint : public T
    {
    public:     // constructor
        using T::T;

    public:     // methods
        void compute(Ref<Vector> rs, const std::vector<Ref<Vector>>& gs, const std::vector<Ref<Matrix>>& hs) const override
        {
            pybind11::gil_scoped_acquire acquire;
            PYBIND11_OVERLOAD_PURE(void, T, compute, rs, gs, hs);
        }
    };

    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Type = Constraint;
        using Trampoline = PyConstraint<Type>;
        using Holder = Pointer<Type>;

        py::class_<Type, Trampoline, Holder>(m, "Constraint")
            // constructors
            .def(py::init<>())
            .def(py::init<index, index>(), "nb_equations"_a, "nb_variables"_a)
            // read-only properties
            .def_property_readonly("nb_equations", &Type::nb_equations)
            .def_property_readonly("nb_variables", &Type::nb_variables)
            // properties
            .def_property("equations", &Type::equations, &Type::set_equations)
            .def_property("is_active", &Type::is_active, &Type::set_active)
            .def_property("name", &Type::name, &Type::set_name)
            .def_property("variables", &Type::variables, &Type::set_variables)
            // methods
            .def("compute", &Type::compute, "fs"_a, "gs"_a, "hs"_a)
            .def("equation", &Type::equation, "index"_a, py::return_value_policy::reference_internal)
            .def("variable", &Type::variable, "index"_a, py::return_value_policy::reference_internal)
        ;
    }
};

} // namespace eqlib