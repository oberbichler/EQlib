#pragma once

#include "Define.h"

namespace EQlib {

class Solver {
public:     // methods
    virtual void analyze_pattern(Ref<const Sparse> a) = 0;

    virtual bool set_matrix(Ref<const Sparse> a) = 0;

    virtual bool solve(Ref<const Vector> b, Ref<Vector> x) = 0;
};

} // namespace EQlib