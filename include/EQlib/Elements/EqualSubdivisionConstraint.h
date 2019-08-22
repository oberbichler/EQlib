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

class EqualSubdivisionConstraint : public Element
{
private:    // types
    using HyperJet = hyperjet::HyperJet<double>;

private:    // variables
    std::vector<Pointer<Node>> m_nodes;
    double m_weight;
    size_t m_nb_dofs;

public:     // constructors
    EqualSubdivisionConstraint(
        const std::vector<Pointer<Node>>& nodes,
        const double weight)
    : m_nodes(nodes)
    , m_weight(weight)
    {
        m_nb_dofs = 3 * nodes.size();
    }

public:     // methods
    std::vector<Pointer<Parameter>> dofs() const override
    {
        std::vector<Pointer<Parameter>> dof_list(m_nb_dofs);

        size_t offset = 0;

        for (const auto& node : m_nodes) {
            dof_list[offset++] = node->x();
            dof_list[offset++] = node->y();
            dof_list[offset++] = node->z();
        }

        assert(offset == m_nb_dofs);

        return dof_list;
    }

    double compute(Ref<Vector> g, Ref<Matrix> h) const override
    {
        using namespace hyperjet;

        HyperJet x_a(m_nodes[0]->x()->act_value(), m_nb_dofs);
        HyperJet y_a(m_nodes[0]->y()->act_value(), m_nb_dofs);
        HyperJet z_a(m_nodes[0]->z()->act_value(), m_nb_dofs);

        size_t offset = 0;

        x_a.g(offset++) = 1;
        y_a.g(offset++) = 1;
        z_a.g(offset++) = 1;

        std::vector<HyperJet> lengths(m_nodes.size() - 1);

        for (size_t i = 1; i < m_nodes.size(); i++) {
            HyperJet x_b(m_nodes[i]->x()->act_value(), m_nb_dofs);
            HyperJet y_b(m_nodes[i]->y()->act_value(), m_nb_dofs);
            HyperJet z_b(m_nodes[i]->z()->act_value(), m_nb_dofs);

            x_b.g(offset++) = 1;
            y_b.g(offset++) = 1;
            z_b.g(offset++) = 1;

            HyperJet dx = x_b - x_a;
            HyperJet dy = y_b - y_a;
            HyperJet dz = z_b - z_a;

            lengths[i - 1] += (dx * dx + dy * dy + dz * dz).sqrt();

            x_a = x_b;
            y_a = y_b;
            z_a = z_b;
        }

        HyperJet total_length(m_nb_dofs);

        for (size_t i = 0; i < lengths.size(); i++) {
            total_length += lengths[i];
        }

        HyperJet average_length = total_length / (m_nodes.size() - 1);

        HyperJet f(m_nb_dofs);

        for (size_t i = 0; i < lengths.size(); i++) {
            f += (lengths[i] - average_length).pow(2);
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

        using Type = EqualSubdivisionConstraint;
        using Holder = Pointer<Type>;
        using Base = Element;

        py::class_<Type, Base, Holder>(m, "EqualSubdivisionConstraint")
            .def(py::init<std::vector<Pointer<Node>>, double>(), "nodes"_a,
            "weight"_a=1.0)
        ;
    }
};

} // namespace EQlib