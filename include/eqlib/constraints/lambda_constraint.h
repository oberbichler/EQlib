#pragma once

#include "../constraint.h"

#include <functional>
#include <vector>

namespace eqlib {
namespace constraints {

struct LambdaConstraint : public Constraint {
    using Fn = std::function<void(Ref<const Vector>, Ref<Vector>, Ref<Matrix>, Ref<Matrix>, const Request)>;

    LambdaConstraint(Equations equations, Parameters parameters, Fn fn)
        : Constraint()
        , _fn(fn)
    {
        set_equations(equations);
        set_parameters(parameters);
    }

    // function

    Fn _fn;

    // compute

    void compute(Ref<Vector> g, Ref<Matrix> dg, Ref<Matrix> hm, const Request request)
    {
        const Vector x = parameter_values();
        _fn(x, g, dg, hm, request);
    }
};

} // namespace constraints
} // namespace eqlib