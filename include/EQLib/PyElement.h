#pragma once

#include "Element.h"

namespace EQLib {

class PyElement
: public Element
{
public:
    using Element::Element;

    std::vector<Dof> dofs() const override
    {
        using ReturnType = std::vector<Dof>;
        PYBIND11_OVERLOAD_PURE(ReturnType, Element, dofs);
    }

    std::pair<Matrix, Vector> compute(py::dict options) const override
    {
        using ReturnType = std::pair<Matrix, Vector>;
        PYBIND11_OVERLOAD_PURE(ReturnType, Element, compute, options);
    }
};

} // namespace EQLib