#pragma once

#include "../Define.h"
#include "../Element.h"
#include "../Node.h"
#include "../Parameter.h"

#include <hyperjet/HyperJet.h>
#include <hyperjet/Jet.h>

#include <variant>
#include <vector>

namespace EQlib {

class BoundaryConstraint : public Element
{
private:    // types
    using HyperJet = hyperjet::HyperJet<double>;

private:    // variables
    Pointer<Parameter> m_parameter;
    double m_lower_bound;
    double m_upper_bound;
    double m_weight;
    size_t m_nb_dofs;

public:     // constructors
    BoundaryConstraint(
        const Pointer<Parameter>& parameter,
        const double lower_bound,
        const double upper_bound,
        const double weight)
    : m_parameter(parameter)
    , m_lower_bound(lower_bound)
    , m_upper_bound(upper_bound)
    , m_weight(weight)
    {
        m_nb_dofs = 1;
    }

public:     // methods
    std::vector<Pointer<Parameter>> dofs() const override
    {
        std::vector<Pointer<Parameter>> dof_list(m_nb_dofs);

        size_t offset = 0;

        dof_list[offset++] = m_parameter;

        assert(offset == m_nb_dofs);

        return dof_list;
    }

    double compute(Ref<Vector> g, Ref<Matrix> h) const override
    {
        using namespace hyperjet;

        HyperJet value(m_parameter->act_value(), m_nb_dofs);

        size_t offset = 0;

        value.g(offset++) = 1;

        HyperJet f(m_nb_dofs);

        if (m_parameter->act_value() < m_lower_bound) {
            f = 0.5 * (value - m_lower_bound).pow(2) * m_weight;
        }

        if (m_parameter->act_value() > m_upper_bound) {
            f = 0.5 * (value - m_upper_bound).pow(2) * m_weight;
        }

        if (g.size() > 0) {
            g = f.g();
        }

        if (h.size() > 0) {
            h = f.h();
        }

        assert(offset == m_nb_dofs);

        return f.f();
    }

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Type = BoundaryConstraint;
        using Holder = Pointer<Type>;
        using Base = Element;

        py::class_<Type, Base, Holder>(m, "BoundaryConstraint")
            .def(py::init<Pointer<Parameter>, double, double, double>(),
            "parameter"_a, "lower_bound"_a, "upper_bound"_a, "weight"_a=1.0)
        ;
    }
};

} // namespace EQlib