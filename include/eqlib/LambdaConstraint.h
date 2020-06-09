#pragma once

#include "Constraint.h"
#include "Define.h"
#include "Equation.h"
#include "Variable.h"

#include <functional>
#include <vector>

namespace eqlib {

class LambdaConstraint : public Constraint {
private: // types
    using Type = LambdaConstraint;
    using ComputeFunction = std::function<void(const std::vector<Pointer<Equation>>&, const std::vector<Pointer<Variable>>&, Ref<Vector>, const std::vector<Ref<Vector>>&, const std::vector<Ref<Matrix>>)>;

private: // variables
    ComputeFunction m_compute;

public: // constructors
    LambdaConstraint(
        const std::vector<Pointer<Equation>>& equations,
        const std::vector<Pointer<Variable>>& variables,
        ComputeFunction compute)
        : Constraint{}
        , m_compute{compute}
    {
        m_equations = equations;
        m_variables = variables;
    }

public: // methods
    void compute(Ref<Vector> rs, const std::vector<Ref<Vector>>& gs, const std::vector<Ref<Matrix>>& hs) const override
    {
        pybind11::gil_scoped_acquire acquire;
        m_compute(m_equations, m_variables, rs, gs, hs);
    }

public: // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Base = Constraint;
        using Holder = Pointer<Type>;

        py::class_<Type, Base, Holder>(m, "LambdaConstraint")
            // constructors
            .def(py::init<const std::vector<Pointer<Equation>>&, const std::vector<Pointer<Variable>>&, ComputeFunction>(), "equations"_a, "variables"_a, "compute"_a);
    }
};

} // namespace eqlib