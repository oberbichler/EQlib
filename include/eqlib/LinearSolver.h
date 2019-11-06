#pragma once

#include "Define.h"

#include <string>

namespace eqlib {

class LinearSolver
{
public:     // methods
    virtual std::string solver_name() const = 0;

    virtual bool analyze(Ref<const Sparse> a) = 0;

    virtual bool factorize(Ref<const Sparse> a) = 0;

    virtual bool solve(Ref<const Sparse> a, Ref<const Vector> b, Ref<Vector> x) = 0;
};

} // namespace eqlib