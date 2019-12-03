#pragma once

#include "Define.h"
#include "Objective.h"
#include "Variable.h"

#include <functional>
#include <vector>

namespace eqlib {

class LambdaObjective : public Objective {
private: // types
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

public: // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Base = Objective;
        using Holder = Pointer<Type>;

        py::class_<Type, Base, Holder>(m, "LambdaObjective")
            // constructors
            .def(py::init<const std::vector<Pointer<Variable>>&, ComputeFunction>(), "variables"_a, "compute"_a);
    }
};

} // namespace eqlib