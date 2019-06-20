#pragma once

#include "Element.h"

namespace EQlib {

class PyElement : public Element
{
public:     // constructor
    using Element::Element;

public:     // getters and setters
    std::vector<Dof> dofs() const override
    {
        using ReturnType = std::vector<Dof>;
        pybind11::gil_scoped_acquire acquire;
        PYBIND11_OVERLOAD_PURE(ReturnType, Element, dofs);
    }

public:     // methods
    double compute(Ref<Vector> g, Ref<Matrix> h) const override
    {
        pybind11::gil_scoped_acquire acquire;
        PYBIND11_OVERLOAD_PURE(double, Element, compute, g, h);
    }

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Type = Element;
        using Trampoline = PyElement;
        using Holder = std::shared_ptr<Type>;

        py::class_<Type, Trampoline, Holder>(m, "Element", py::dynamic_attr())
            .def(py::init<>())
            .def("dofs", &Type::dofs)
            .def("compute", &Type::compute, "g"_a, "h"_a)
        ;
    }
};

} // namespace EQlib