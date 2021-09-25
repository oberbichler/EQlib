#pragma once

#include "common.h"
#include "objective.h"
#include "variable.h"

#include <vector>
#include <functional>

namespace eqlib
{

    struct LambdaObjective : public Objective
    {
        using Fn = std::function<double(Ref<const Vector>, Ref<Vector>, Ref<Matrix>, const Request)>;

        LambdaObjective(Variables variables, Fn fn) : Objective(), _fn(fn)
        {
            set_variables(variables);
        }

        // function

        Fn _fn;

        // compute

        double compute(Ref<Vector> df, Ref<Matrix> hm, const Request request)
        {
            const Vector x = variable_values();
            return _fn(x, df, hm, request);
        }
    };

} // namespace eqlib