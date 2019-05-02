#pragma once

#include <EQlib/Element.h>
#include <EQlib/Node.h>

#include <HyperJet/HyperJet.h>

#include <vector>

namespace EQlib {

class LocationConstraint : public Element
{
private:    // variables
    std::vector<std::shared_ptr<Node>> m_nodes;
    Matrix m_shape_functions;
    Vector3D m_target;
    double m_penalty;

public:     // constructor
    LocationConstraint(
        std::vector<std::shared_ptr<Node>> nodes,
        Matrix shape_functions,
        Vector3D target,
        double penalty)
    : m_nodes(nodes)
    , m_shape_functions(shape_functions)
    , m_target(target)
    , m_penalty(penalty)
    { }

public:     // methods
    Eigen::Matrix<HyperJet::HyperJet<double>, 3, 1> act_evaluate() const
    {
        const size_t nb_dofs = m_nodes.size() * 3;

        Vector3D xyz = Vector3D::Zero();
        Vector dx = Vector::Zero(nb_dofs);
        Vector dy = Vector::Zero(nb_dofs);
        Vector dz = Vector::Zero(nb_dofs);

        for (size_t i = 0; i < m_nodes.size(); i++) {
            xyz += m_nodes[i]->act_location() * m_shape_functions(0, i);
            dx[i * 3 + 0] = m_shape_functions(0, i);
            dy[i * 3 + 1] = m_shape_functions(0, i);
            dz[i * 3 + 2] = m_shape_functions(0, i);
        }

        Eigen::Matrix<HyperJet::HyperJet<double>, 3, 1> result;

        result[0] = HyperJet::HyperJet<double>(xyz(0), dx);
        result[1] = HyperJet::HyperJet<double>(xyz(1), dy);
        result[2] = HyperJet::HyperJet<double>(xyz(2), dz);

        return result;
    }

    std::vector<Dof> dofs() const override
    {
        std::vector<Dof> dof_list(m_nodes.size() * 3);

        for (size_t i = 0; i < m_nodes.size(); i++) {
            dof_list[i * 3 + 0] = m_nodes[i]->x().dof();
            dof_list[i * 3 + 1] = m_nodes[i]->y().dof();
            dof_list[i * 3 + 2] = m_nodes[i]->z().dof();
        }

        return dof_list;
    }

    std::tuple<double, Vector, Matrix> compute() const override
    {
        const auto location = act_evaluate();

        const auto v = m_target - location;

        const auto f = 0.5 * v.dot(v);

        return {f.f(), f.g(), f.h()};
    }
};

} // namespace EQlib