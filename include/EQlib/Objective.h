#pragma once

#include "Define.h"
#include "Variable.h"

#include <vector>

namespace EQlib {

class Objective
{
private:    // variables
    bool m_is_active;

public:     // methods
    Objective() : m_is_active(true) { }

    virtual std::vector<Pointer<Variable>> variables() const = 0;

    virtual double compute(Ref<Vector> g, Ref<Matrix> h) const = 0;

    bool is_active() const noexcept
    {
        return m_is_active;
    }

    void set_is_active(const bool value) noexcept
    {
        m_is_active = value;
    }

public:     // python
    template <typename T>
    class PyObjective : public T
    {
    public:     // constructor
        using T::T;

    public:     // getters and setters
        std::vector<Pointer<Variable>> variables() const override
        {
            using ReturnType = std::vector<Pointer<Variable>>;
            pybind11::gil_scoped_acquire acquire;
            PYBIND11_OVERLOAD_PURE(ReturnType, T, variables);
        }

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
            .def(py::init<>())
            .def("variables", &Type::variables)
            .def("compute", &Type::compute, "g"_a, "h"_a)
            .def_property("is_active", &Type::is_active, &Type::set_is_active)
        ;
    }
};

} // namespace EQlib