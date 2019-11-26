#pragma once

#include "Define.h"
#include "Variable.h"

#include <string>
#include <vector>

namespace eqlib {

class Objective
{
protected:  // variables
    std::string m_name;
    bool m_is_active;
    std::vector<Pointer<Variable>> m_variables;

public:     // constructors
    Objective() : m_is_active(true), m_name("")
    {
    }

    Objective(const index nb_variables)
    : m_variables(nb_variables), m_is_active(true), m_name("")
    {
    }

    virtual ~Objective() = default;

public:     // methods
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

protected:  // methods
    void set_variables(const std::vector<Pointer<Variable>>& value)
    {
        m_variables = value;
    }

public:     // python
    template <typename T>
    class PyObjective : public T
    {
    public:     // constructor
        using T::T;

    public:     // methods
        double compute(Ref<Vector> g, Ref<Matrix> h) const override
        {
            pybind11::gil_scoped_acquire acquire;
            PYBIND11_OVERLOAD_PURE(double, T, compute, g, h);
        }
    };

    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Type = Objective;
        using Trampoline = PyObjective<Type>;
        using Holder = Pointer<Type>;

        py::class_<Type, Trampoline, Holder>(m, "Objective")
            // constructors
            .def(py::init<>())
            .def(py::init<index>(), "nb_variables"_a)
            // read-only properties
            .def_property_readonly("nb_variables", &Type::nb_variables)
            // properties
            .def_property("is_active", &Type::is_active, &Type::set_active)
            .def_property("name", &Type::name, &Type::set_name)
            .def_property("variables", &Type::variables, &Type::set_variables)
            // methods
            .def("compute", &Type::compute, "g"_a, "h"_a)
        ;
    }
};

} // namespace eqlib