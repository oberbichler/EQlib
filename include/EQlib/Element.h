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

    virtual std::pair<Matrix, Vector> compute() const = 0;
};

} // namespace EQlib