#pragma once

#include "Define.h"
#include "Objective.h"
#include "Variable.h"

#include <functional>
#include <vector>

namespace eqlib {

class LambdaObjective : public Objective {
public: // types
    using Type = LambdaObjective;
    using ComputeFunction = std::function<double(const std::vector<Pointer<Variable>>&, Ref<Vector>, Ref<Matrix>)>;

private: // variables
    ComputeFunction m_compute;

public: // constructors
    LambdaObjective(
        const std::vector<Pointer<Variable>>& variables,
        ComputeFunction compute)
        : Objective{}
        , m_compute{compute}
    {
        m_variables = variables;
    }

public: // methods
    double compute(Ref<Vector> g, Ref<Matrix> h) const override
    {
        return m_compute(m_variables, g, h);
    }
};

} // namespace eqlib
