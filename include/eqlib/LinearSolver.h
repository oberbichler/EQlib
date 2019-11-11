#pragma once

#include "Define.h"

#include <string>

namespace eqlib {

class LinearSolver
{
public:     // methods
    virtual std::string solver_name() const = 0;

    virtual bool analyze(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a) = 0;

    virtual bool factorize(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a) = 0;

    virtual bool solve(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a, Ref<const Vector> b, Ref<Vector> x) = 0;
};

} // namespace eqlib