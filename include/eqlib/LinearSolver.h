#pragma once

#include "Define.h"

#include <string>

namespace eqlib {

class LinearSolver {
private: // types
    using Type = LinearSolver;

private: // variables
    std::string m_solver_name;

public: // constructors
    virtual ~LinearSolver() = default;

public: // methods
    std::string solver_name() const
    {
        return m_solver_name;
    }

    void set_solver_name(const std::string& value)
    {
        m_solver_name = value;
    }

    virtual bool analyze(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a)
    {
        return false;
    }

    virtual bool factorize(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a)
    {
        return false;
    }

    virtual bool solve(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a, Ref<const Vector> b, Ref<Vector> x) = 0;
};

} // namespace eqlib
