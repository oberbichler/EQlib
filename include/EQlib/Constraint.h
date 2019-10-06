#pragma once

#include "Define.h"
#include "Equation.h"
#include "Variable.h"

#include <vector>

namespace EQlib {

class Constraint
{
private:    // variables
    bool m_is_active;
    std::vector<Pointer<Equation>> m_equations;
    std::vector<Pointer<Variable>> m_variables;

public:     // methods
    Constraint() : m_is_active(true) { }

    virtual ~Constraint() = default;

    const Pointer<Equation>& equation(const index i) const
    {
        return m_equations[i];
    }

    const Pointer<Variable>& variable(const index i) const
    {
        return m_variables[i];
    }

    const std::vector<Pointer<Equation>>& equations() const
    {
        return m_equations;
    }

    const std::vector<Pointer<Variable>>& variables() const
    {
        return m_variables;
    }

    virtual void compute(Ref<Vector> rs, const std::vector<Ref<Vector>>& gs,
        const std::vector<Ref<Matrix>>& hs) const = 0;

    bool is_active() const noexcept
    {
        return m_is_active;
    }

    void set_active(const bool value) noexcept
    {
        m_is_active = value;
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
        void compute(Ref<Vector> rs, const std::vector<Ref<Vector>>& gs,
            const std::vector<Ref<Matrix>>& hs) const override
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
            .def(py::init<>())
            .def_property("equations", &Type::equations, &Type::set_equations)
            .def_property("variables", &Type::variables, &Type::set_variables)
            .def("compute", &Type::compute, "fs"_a, "gs"_a, "hs"_a)
            .def_property("is_active", &Type::is_active, &Type::set_active)
        ;
    }
};

} // namespace EQlib