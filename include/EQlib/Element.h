#pragma once

#include "Define.h"
#include "Dof.h"

#include <vector>

namespace EQlib {

class Element
{
public:     // methods
    Element() { }

    virtual std::vector<Dof> dofs() const = 0;

    virtual double compute(Ref<Vector> g, Ref<Matrix> h) const = 0;
};

} // namespace EQlib