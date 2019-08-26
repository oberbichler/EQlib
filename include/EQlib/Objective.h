#pragma once

#include "Define.h"
#include "Variable.h"

#include <vector>

namespace EQlib {

class Objective
{
public:     // methods
    Objective() { }

    virtual std::vector<Pointer<Variable>> variables() const = 0;

    virtual double compute(Ref<Vector> g, Ref<Matrix> h) const = 0;

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
        ;
    }
};

} // namespace EQlib