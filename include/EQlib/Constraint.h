#pragma once

#include "Define.h"
#include "Equation.h"
#include "Variable.h"

#include <vector>

namespace EQlib {

class Constraint
{
public:     // methods
    Constraint() { }

    virtual std::vector<Pointer<Equation>> equations() const = 0;

    virtual std::vector<Pointer<Variable>> variables() const = 0;

    virtual void compute(Ref<Vector> rs, const std::vector<Ref<Vector>>& gs,
        const std::vector<Ref<Matrix>>& hs) const = 0;

public:     // python
    template <typename T>
    class PyConstraint : public T
    {
    public:     // constructor
        using T::T;

    public:     // getters and setters
        std::vector<Pointer<Equation>> equations() const override
        {
            using ReturnType = std::vector<Pointer<Equation>>;
            pybind11::gil_scoped_acquire acquire;
            PYBIND11_OVERLOAD_PURE(ReturnType, T, equations);
        }

        std::vector<Pointer<Variable>> variables() const override
        {
            using ReturnType = std::vector<Pointer<Variable>>;
            pybind11::gil_scoped_acquire acquire;
            PYBIND11_OVERLOAD_PURE(ReturnType, T, variables);
        }

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
            .def("equations", &Type::equations)
            .def("variables", &Type::variables)
            .def("compute", &Type::compute, "fs"_a, "gs"_a, "hs"_a)
        ;
    }
};

} // namespace EQlib