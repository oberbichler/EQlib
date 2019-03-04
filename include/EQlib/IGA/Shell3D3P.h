#pragma once

#include <Eigen/Geometry>

#include <EQlib/Element.h>
#include <EQlib/Node.h>

#include <HyperJet/HyperJet.h>

#include <vector>

namespace EQlib {

class Shell3D3P : public Element
{
    using HyperDual = HyperJet::HyperJet<double>;

    template <int TSize>
    using HyperDualVector = Eigen::Matrix<HyperDual, TSize, 1>;

    template <int TRows, int TCols>
    using HyperDualMatrix = Eigen::Matrix<HyperDual, TRows, TCols>;

private:    // variables
    std::vector<std::shared_ptr<Node>> m_nodes;
    Matrix m_shape_functions;
    double m_thickness;
    double m_young_modulus;
    double m_poisson_ratio;
    double m_weight;

public:     // constructor
    Shell3D3P(
        std::vector<std::shared_ptr<Node>> nodes,
        Matrix shape_functions,
        double thickness,
        double young_modulus,
        double poisson_ratio,
        double weight)
    : m_nodes(nodes)
    , m_shape_functions(shape_functions)
    , m_thickness(thickness)
    , m_young_modulus(young_modulus)
    , m_poisson_ratio(poisson_ratio)
    , m_weight(weight)
    { }

    Vector3D ref_evaluate(const int index) const
    {
        const size_t nb_dofs = m_nodes.size() * 3;

        Vector3D xyz = Vector3D::Zero();

        for (size_t i = 0; i < m_nodes.size(); i++) {
            xyz += m_nodes[i]->ref_location() * m_shape_functions(index, i);
        }

        return xyz;
    }

    HyperDualVector<3> act_evaluate(const int index) const
    {
        const size_t nb_dofs = m_nodes.size() * 3;

        Vector3D xyz = Vector3D::Zero();

        HyperDualVector<3> result;

        result(0) = HyperDual(xyz(0), nb_dofs);
        result(1) = HyperDual(xyz(1), nb_dofs);
        result(2) = HyperDual(xyz(2), nb_dofs);

        for (size_t i = 0; i < m_nodes.size(); i++) {
            const double shape_function = m_shape_functions(index, i);
            xyz += m_nodes[i]->act_location() * shape_function;
            result[0].g()(i * 3 + 0) = shape_function;
            result[1].g()(i * 3 + 1) = shape_function;
            result[2].g()(i * 3 + 2) = shape_function;
        }

        result(0).f() = xyz(0);
        result(1).f() = xyz(1);
        result(2).f() = xyz(2);

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

    std::pair<Matrix, Vector> compute(py::dict options) const override
    {
        Eigen::Matrix3d Dm;
        Dm << 1.0, m_poisson_ratio, 0,
              m_poisson_ratio, 1.0, 0,
              0, 0, (1.0 - m_poisson_ratio) / 2.0;
        Dm *= m_young_modulus * m_thickness / (1.0 - std::pow(m_poisson_ratio, 2));

        Eigen::Matrix3d Db;
        Db << 1.0, m_poisson_ratio, 0,
              m_poisson_ratio, 1.0, 0,
              0, 0, (1.0 - m_poisson_ratio) / 2.0;
        Db *= m_young_modulus * std::pow(m_thickness, 3) / (12.0 * (1.0 - std::pow(m_poisson_ratio, 2)));

        // reference configuration

        const auto A1 = ref_evaluate(1);
        const auto A2 = ref_evaluate(2);

        const auto A1_1 = ref_evaluate(3);
        const auto A1_2 = ref_evaluate(4);
        const auto A2_2 = ref_evaluate(5);

        const auto A11 = A1.dot(A1);
        const auto A22 = A2.dot(A2);
        const auto A12 = A1.dot(A2);

        auto A3 = A1.cross(A2);
        const auto dA = A3.norm();
        A3 = A3 / dA;

        const auto B11 = A1_1.dot(A3);
        const auto B12 = A1_2.dot(A3);
        const auto B22 = A2_2.dot(A3);

        // transformation

        Vector3D e1 = A1 / A1.norm();
        Vector3D e2 = A2 - A2.dot(e1) * e1;
        e2 /= e2.norm();

        double det = A11 * A22 - A12 * A12;

        Vector3D g_ab_con;
        g_ab_con << A22 / det, A11 / det, -A12 / det;

        Vector3D g_con1 = g_ab_con[0] * A1 + g_ab_con[2] * A2;
        Vector3D g_con2 = g_ab_con[2] * A1 + g_ab_con[1] * A2;

        const double eg11 = e1.dot(g_con1);
        const double eg12 = e1.dot(g_con2);
        const double eg21 = e2.dot(g_con1);
        const double eg22 = e2.dot(g_con2);

        Eigen::Matrix3d Tm;
        Tm << eg11 * eg11, eg12 * eg12, 2 * eg11 * eg12,
              eg21 * eg21, eg22 * eg22, 2 * eg21 * eg22,
              2 * eg11 * eg21, 2 * eg12 * eg22, 2 * (eg11 * eg22 + eg12 * eg21);

        // actual configuration

        const HyperDualVector<3> a1 = act_evaluate(1);
        const HyperDualVector<3> a2 = act_evaluate(2);

        const HyperDualVector<3> a1_1 = act_evaluate(3);
        const HyperDualVector<3> a1_2 = act_evaluate(4);
        const HyperDualVector<3> a2_2 = act_evaluate(5);

        const HyperDual a11 = a1.dot(a1);
        const HyperDual a22 = a2.dot(a2);
        const HyperDual a12 = a1.dot(a2);

        HyperDualVector<3> a3 = a1.cross(a2);
        a3 = a3 / a3.norm();

        const HyperDual b11 = a1_1.dot(a3);
        const HyperDual b12 = a1_2.dot(a3);
        const HyperDual b22 = a2_2.dot(a3);

        const HyperDualVector<3> eps = Tm * 0.5 * HyperDualVector<3>(a11 - A11, a22 - A22, a12 - A12);
        const HyperDualVector<3> kap = Tm *       HyperDualVector<3>(B11 - b11, B12 - b12, B22 - b22);

        const HyperDualMatrix<1, 1> f = (0.5 * (eps.transpose() * Dm * eps +
                                                kap.transpose() * Db * kap));

        return {f(0, 0).h(), -f(0, 0).g()};
    }
};

} // namespace EQlib