#pragma once

#include "common.h"

#include <string>

namespace eqlib {

struct LinearSolver {
    virtual ~LinearSolver() = default;

    // solver_name

    std::string m_solver_name;

    std::string solver_name() const
    {
        return m_solver_name;
    }

    void set_solver_name(const std::string& value)
    {
        m_solver_name = value;
    }

    // analyze

    virtual bool analyze(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a)
    {
        return false;
    }

    // factorize

    virtual bool factorize(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a)
    {
        return false;
    }

    // solve

    virtual bool solve(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a, Ref<const Vector> b, Ref<Vector> x) = 0;
};

} // namespace eqlib