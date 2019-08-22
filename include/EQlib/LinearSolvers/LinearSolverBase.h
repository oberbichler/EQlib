#pragma once

#include "../Define.h"

namespace EQlib {

struct LinearSolverBase
{
    virtual void analyze(Ref<const Sparse> a) = 0;

    virtual void factorize(Ref<const Sparse> a) = 0;

    virtual void solve(Ref<const Vector> b, Ref<Vector> x) = 0;

    virtual Eigen::ComputationInfo info() const = 0;
};

} // namespace EQlib