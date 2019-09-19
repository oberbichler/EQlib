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

class EqualityConstraint : public Element
{
private:    // types
    using HyperJet = hyperjet::HyperJet<double>;

private:    // variables
    std::vector<Pointer<Parameter>> m_parameters;
    double m_weight;
    size_t m_nb_dofs;

public:     // constructors
    EqualityConstraint(
        const std::vector<Pointer<Parameter>>& parameters,
        const double weight)
    : m_parameters(parameters)
    , m_weight(weight)
    {
        m_nb_dofs = parameters.size();
    }

public:     // methods
    std::vector<Pointer<Parameter>> dofs() const override
    {
        std::vector<Pointer<Parameter>> dof_list(m_nb_dofs);

        size_t offset = 0;

        for (const auto& parameter : m_parameters) {
            dof_list[offset++] = parameter;
        }

        assert(offset == m_nb_dofs);

        return dof_list;
    }

    double compute(Ref<Vector> g, Ref<Matrix> h) const override
    {
        using namespace hyperjet;

        size_t offset = 0;

        std::vector<HyperJet> values(m_parameters.size());

        for (size_t i = 0; i < m_parameters.size(); i++) {
            HyperJet value(m_parameters[i]->act_value(), m_nb_dofs);

            value.g(offset++) = 1;

            values[i] = value;
        }

        HyperJet sum(m_nb_dofs);

        for (size_t i = 0; i < values.size(); i++) {
            sum += values[i];
        }

        HyperJet average = sum / m_parameters.size();

        HyperJet f(m_nb_dofs);

        for (size_t i = 0; i < values.size(); i++) {
            f += (values[i] - average).pow(2);
        }

        f = 0.5 * f * m_weight;

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

        using Type = EqualityConstraint;
        using Holder = Pointer<Type>;
        using Base = Element;

        py::class_<Type, Base, Holder>(m, "EqualityConstraint")
            .def(py::init<std::vector<Pointer<Parameter>>, double>(),
                "parameters"_a, "weight"_a=1.0)
        ;
    }
};

} // namespace EQlib