#pragma once

#include "Define.h"
#include "Parameter.h"

#include <vector>

namespace EQlib {

class Element
{
public:     // methods
    Element() { }

    virtual ~Element() = default;

    virtual std::vector<Pointer<Parameter>> dofs() const = 0;

    virtual double compute(Ref<Vector> g, Ref<Matrix> h) const = 0;

public:     // python
    template <typename T>
    class PyElement : public T
    {
    public:     // constructor
        using T::T;

    public:     // getters and setters
        std::vector<Pointer<Parameter>> dofs() const override
        {
            using ReturnType = std::vector<Pointer<Parameter>>;
            pybind11::gil_scoped_acquire acquire;
            PYBIND11_OVERLOAD_PURE(ReturnType, T, dofs);
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

        using Type = Element;
        using Trampoline = PyElement<Type>;
        using Holder = Pointer<Type>;

        py::class_<Type, Trampoline, Holder>(m, "Element", py::dynamic_attr())
            .def(py::init<>())
            .def("dofs", &Type::dofs)
            .def("compute", &Type::compute, "g"_a, "h"_a)
            .def("compute", [](const Type& self) {
                Vector g(0);
                Matrix h(0, 0);
                return self.compute(g, h);
            })
        ;
    }
};

} // namespace EQlib