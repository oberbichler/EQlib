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
    std::vector<Dof> dofs() const override
    {
        std::vector<Dof> dof_list(m_nb_dofs);

        size_t offset = 0;

        dof_list[offset++] = m_node->x();
        dof_list[offset++] = m_node->y();
        dof_list[offset++] = m_node->z();

        for (const auto& [force, node] : m_connections) {
            dof_list[offset++] = node->x();
            dof_list[offset++] = node->y();
            dof_list[offset++] = node->z();
            dof_list[offset++] = *force;
        }

        for (const auto& [force, direction] : m_loads) {
            if (std::holds_alternative<Pointer<Parameter>>(force)) {
                dof_list[offset++] = *std::get<Pointer<Parameter>>(force);
            }
        }

        assert(offset == m_nb_dofs);

        return dof_list;
    }

    double compute(Ref<Vector> g, Ref<Matrix> h) const override
    {
        using namespace hyperjet;
        
        HyperJet x(m_node->x(), m_nb_dofs);
        HyperJet y(m_node->y(), m_nb_dofs);
        HyperJet z(m_node->z(), m_nb_dofs);

        size_t offset = 0;

        x.g(offset++) = 1;
        y.g(offset++) = 1;
        z.g(offset++) = 1;

        Eigen::Matrix<HyperJet, 3, 1> residual;

        residual[0] = HyperJet(m_nb_dofs);
        residual[1] = HyperJet(m_nb_dofs);
        residual[2] = HyperJet(m_nb_dofs);

        for (const auto& [force, node] : m_connections) {
            Eigen::Matrix<HyperJet, 3, 1> f_i;

            f_i[0] = HyperJet(node->x(), m_nb_dofs);
            f_i[1] = HyperJet(node->y(), m_nb_dofs);
            f_i[2] = HyperJet(node->z(), m_nb_dofs);

            HyperJet s_i(*force, m_nb_dofs);

            f_i[0].g(offset++) = 1;
            f_i[1].g(offset++) = 1;
            f_i[2].g(offset++) = 1;

            f_i[0] -= x;
            f_i[1] -= y;
            f_i[2] -= z;

            s_i.g(offset++) = 1;

            const HyperJet l_i = f_i.norm();

            f_i *= s_i / l_i;

            residual += f_i;
        }

        for (const auto& [force, direction] : m_loads) {
            if (std::holds_alternative<Pointer<Parameter>>(force)) {
                HyperJet s_i(*std::get<Pointer<Parameter>>(force), m_nb_dofs);

                s_i.g(offset++) = 1;

                residual += s_i * direction / direction.norm();
            } else {
                const double s_i = std::get<double>(force);

                const Vector3D f_i = s_i * direction / direction.norm();

                residual = residual + f_i;
            }
        }

        const auto f = 0.5 * residual.dot(residual);

        if (h.size() > 0) {
            g = f.g();
            h = f.h();
        } else if (g.size() > 0) {
            g = f.g();
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