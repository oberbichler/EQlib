#pragma once

#include "Define.h"
#include "Dof.h"

#include <vector>

namespace EQLib {

class Element
{
public:     // methods
    Element() { }

    virtual std::vector<Dof> dofs() const = 0;

    virtual std::pair<Matrix, Vector> compute(py::dict options) const = 0;
};

} // namespace EQlib