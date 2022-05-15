#pragma once

#include "../objective.h"

#include <functional>
#include <vector>

namespace eqlib {
namespace objectives {

struct LambdaObjective : public Objective {
    using Fn = std::function<double(Ref<const Vector>, Ref<Vector>, Ref<Matrix>, const Request)>;

    LambdaObjective(Parameters parameters, Fn fn)
        : Objective()
        , _fn(fn)
    {
        set_parameters(parameters);
    }

    // function

    Fn _fn;

    // compute

    double compute(Ref<Vector> df, Ref<Matrix> hm, const Request request)
    {
        const Vector x = parameter_values();
        return _fn(x, df, hm, request);
    }
};

} // namespace objectives
} // namespace eqlib