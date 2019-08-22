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

class NodalEquilibrium : public Element
{
private:    // types
    using HyperJet = hyperjet::HyperJet<double>;
    using Connection = std::pair<Pointer<Parameter>, Pointer<Node>>;
    using Load = std::pair<std::variant<double, Pointer<Parameter>>, Vector3D>;

private:    // variables
    Pointer<Node> m_node;
    std::vector<Connection> m_connections;
    std::vector<Load> m_loads;
    double m_weight;
    size_t m_nb_dofs;

public:     // constructors
    NodalEquilibrium(
        const Pointer<Node>& node,
        const std::vector<Connection>& connections,
        const std::vector<Load>& loads,
        const double weight)
    : m_node(node)
    , m_connections(connections)
    , m_loads(loads)
    , m_weight(weight)
    {
        m_nb_dofs = 3 + 4 * connections.size();

        for (const auto& [force, direction] : m_loads) {
            if (std::holds_alternative<Pointer<Parameter>>(force)) {
                m_nb_dofs++;
            }
        }
    }

public:     // methods
    std::vector<Pointer<Parameter>> dofs() const override
    {
        std::vector<Pointer<Parameter>> dof_list(m_nb_dofs);

        size_t dof = 0;

        dof_list[dof++] = m_node->x();
        dof_list[dof++] = m_node->y();
        dof_list[dof++] = m_node->z();

        for (const auto& [force, node] : m_connections) {
            dof_list[dof++] = node->x();
            dof_list[dof++] = node->y();
            dof_list[dof++] = node->z();
            dof_list[dof++] = *force;
        }

        for (const auto& [force, direction] : m_loads) {
            if (std::holds_alternative<Pointer<Parameter>>(force)) {
                dof_list[dof++] = *std::get<Pointer<Parameter>>(force);
            }
        }

        assert(dof == m_nb_dofs);

        return dof_list;
    }

    double compute(Ref<Vector> g, Ref<Matrix> h) const override
    {
        using namespace hyperjet;

        size_t dof = 0;

        const auto x = HyperJet::variable(m_node->x(), m_nb_dofs, dof++);
        const auto y = HyperJet::variable(m_node->y(), m_nb_dofs, dof++);
        const auto z = HyperJet::variable(m_node->z(), m_nb_dofs, dof++);

        Eigen::Matrix<HyperJet, 3, 1> residual;

        residual[0] = HyperJet(m_nb_dofs);
        residual[1] = HyperJet(m_nb_dofs);
        residual[2] = HyperJet(m_nb_dofs);

        for (const auto& [force, node] : m_connections) {
            Eigen::Matrix<HyperJet, 3, 1> direction;

            direction[0] = HyperJet::variable(node->x(), m_nb_dofs, dof++) - x;
            direction[1] = HyperJet::variable(node->y(), m_nb_dofs, dof++) - y;
            direction[2] = HyperJet::variable(node->z(), m_nb_dofs, dof++) - z;

            const auto s = HyperJet::variable(*force, m_nb_dofs, dof++);

            residual = residual + s * direction / direction.norm();
        }

        for (const auto& [force, direction] : m_loads) {
            if (std::holds_alternative<Pointer<Parameter>>(force)) {
                const auto s = HyperJet::variable(
                    *std::get<Pointer<Parameter>>(force), m_nb_dofs, dof++);

                residual = residual + s * direction / direction.norm();
            } else {
                const auto s = std::get<double>(force);

                residual = residual + s * direction / direction.norm();;
            }
        }

        const auto f = 0.5 * residual.dot(residual);

        if (g.size() > 0) {
            g = f.g();
        }

        if (h.size() > 0) {
            h = f.h();
        }

        assert(dof == m_nb_dofs);

        return f.f();
    }

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Type = NodalEquilibrium;
        using Base = Element;
        using Holder = Pointer<Type>;

        py::class_<Type, Base, Holder>(m, "NodalEquilibrium")
            .def(py::init<Pointer<Node>, std::vector<Connection>,
            std::vector<Load>, double>(), "node"_a, "connections"_a, "loads"_a,
            "weight"_a=1.0)
        ;
    }
};

} // namespace EQlib