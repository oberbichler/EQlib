#pragma once

#include <Eigen/Geometry>

#include <EQlib/Objective.h>
#include <EQlib/Point.h>
#include <EQlib/Variable.h>

#include <hyperjet/hyperjet.h>

#include <vector>

namespace EQlib {

class IgaShell3PAD : public Objective
{
private:    // types
    using Type = IgaShell3PAD;

    using Jet = hyperjet::Jet<double>;
    using HyperJet = hyperjet::HyperJet<double>;

    using Jet3D = Eigen::Matrix<Jet, 3, 1>;
    using HyperJet3D = Eigen::Matrix<HyperJet, 3, 1>;

private:    // variables
    std::vector<Pointer<Point>> m_nodes;
    Matrix m_shape_functions;
    double m_thickness;
    double m_young_modulus;
    double m_poisson_ratio;
    Eigen::Matrix3d m_dm;
    Eigen::Matrix3d m_db;
    double ref_a11;
    double ref_a12;
    double ref_a22;
    double ref_b11;
    double ref_b12;
    double ref_b22;
    Eigen::Matrix3d tm;
    double m_weight;
    double m_da;

    template <int TOrder>
    auto ref_geometry(index i) const
    {
        Vector3D value = Vector3D::Zero();

        for (index j = 0; j < length(m_nodes); j++) {
            value += m_nodes[j]->ref_location() * m_shape_functions(i, j);
        }

        if constexpr(TOrder == 0) {
            return value;
        }

        if constexpr(TOrder == 1) {
            Jet3D jet;

            for (index k = 0; k < 3; k++) {
                jet[k] = Jet(value(k), length(m_nodes) * 3);
            }

            for (index j = 0; j < length(m_nodes); j++) {
                jet(0).g(j * 3 + 0) = m_shape_functions(i, j);
                jet(1).g(j * 3 + 1) = m_shape_functions(i, j);
                jet(2).g(j * 3 + 2) = m_shape_functions(i, j);
            }

            return jet;
        }

        if constexpr(TOrder == 2) {
            HyperJet3D hyper_jet;

            for (index k = 0; k < 3; k++) {
                hyper_jet[k] = HyperJet(value(k), length(m_nodes) * 3);
            }

            for (index j = 0; j < length(m_nodes); j++) {
                hyper_jet(0).g(j * 3 + 0) = m_shape_functions(i, j);
                hyper_jet(1).g(j * 3 + 1) = m_shape_functions(i, j);
                hyper_jet(2).g(j * 3 + 2) = m_shape_functions(i, j);
            }

            return hyper_jet;
        }
    }

    template <int TOrder>
    auto act_geometry(index i) const
    {
        Vector3D value = Vector3D::Zero();

        for (index j = 0; j < length(m_nodes); j++) {
            value += m_nodes[j]->act_location() * m_shape_functions(i, j);
        }

        if constexpr(TOrder == 0) {
            return value;
        }

        if constexpr(TOrder == 1) {
            Jet3D jet;

            for (index k = 0; k < 3; k++) {
                jet[k] = Jet(value(k), length(m_nodes) * 3);
            }

            for (index j = 0; j < length(m_nodes); j++) {
                jet(0).g(j * 3 + 0) = m_shape_functions(i, j);
                jet(1).g(j * 3 + 1) = m_shape_functions(i, j);
                jet(2).g(j * 3 + 2) = m_shape_functions(i, j);
            }

            return jet;
        }

        if constexpr(TOrder == 2) {
            HyperJet3D hyper_jet;

            for (index k = 0; k < 3; k++) {
                hyper_jet[k] = HyperJet(value(k), length(m_nodes) * 3);
            }

            for (index j = 0; j < length(m_nodes); j++) {
                hyper_jet(0).g(j * 3 + 0) = m_shape_functions(i, j);
                hyper_jet(1).g(j * 3 + 1) = m_shape_functions(i, j);
                hyper_jet(2).g(j * 3 + 2) = m_shape_functions(i, j);
            }

            return hyper_jet;
        }
    }

public:     // constructor
    IgaShell3PAD(
        std::vector<Pointer<Point>> nodes,
        Matrix shape_functions,
        double thickness,
        double young_modulus,
        double poisson_ratio,
        double weight)
    : m_shape_functions(shape_functions)
    , m_nodes(nodes)
    , m_weight(weight)
    {
        m_dm << 1.0, poisson_ratio, 0,
                poisson_ratio, 1.0, 0,
                0, 0, (1.0 - poisson_ratio) / 2.0;
        m_dm *= young_modulus * thickness / (1.0 - std::pow(poisson_ratio, 2));

        m_db << 1.0, poisson_ratio, 0,
                poisson_ratio, 1.0, 0,
                0, 0, (1.0 - poisson_ratio) / 2.0;
        m_db *= young_modulus * std::pow(thickness, 3) / (12.0 * (1.0 - std::pow(poisson_ratio, 2)));

        m_variables.resize(length(nodes) * 3);
        for (index i = 0; i < length(nodes); i++) {
            m_variables[i * 3 + 0] = nodes[i]->x();
            m_variables[i * 3 + 1] = nodes[i]->y();
            m_variables[i * 3 + 2] = nodes[i]->z();
        }

        // reference configuration

        const Vector3D ref_a1 = ref_geometry<0>(1);
        const Vector3D ref_a2 = ref_geometry<0>(2);

        const Vector3D ref_a1_1 = ref_geometry<0>(3);
        const Vector3D ref_a1_2 = ref_geometry<0>(4);
        const Vector3D ref_a2_2 = ref_geometry<0>(5);

        Vector3D ref_a3 = ref_a1.cross(ref_a2);
        const double ref_da = ref_a3.norm();
        ref_a3 = ref_a3 / ref_da;

        ref_a11 = ref_a1.dot(ref_a1);
        ref_a22 = ref_a2.dot(ref_a2);
        ref_a12 = ref_a1.dot(ref_a2);

        ref_b11 = ref_a1_1.dot(ref_a3);
        ref_b12 = ref_a1_2.dot(ref_a3);
        ref_b22 = ref_a2_2.dot(ref_a3);

        const Vector3D e1 = ref_a1 / ref_a1.norm();
        Vector3D e2 = ref_a2 - ref_a2.dot(e1) * e1;
        e2 /= e2.norm();

        const double det = ref_a11 * ref_a22 - ref_a12 * ref_a12;

        Vector3D g_ab_con;
        g_ab_con << ref_a22 / det, ref_a11 / det, -ref_a12 / det;

        Vector3D g_con1 = g_ab_con[0] * ref_a1 + g_ab_con[2] * ref_a2;
        Vector3D g_con2 = g_ab_con[2] * ref_a1 + g_ab_con[1] * ref_a2;

        const double eg11 = e1.dot(g_con1);
        const double eg12 = e1.dot(g_con2);
        const double eg21 = e2.dot(g_con1);
        const double eg22 = e2.dot(g_con2);

        // Eigen::Matrix3d tm;
        tm << eg11 * eg11, eg12 * eg12, 2 * eg11 * eg12,
              eg21 * eg21, eg22 * eg22, 2 * eg21 * eg22,
              2 * eg11 * eg21, 2 * eg12 * eg22, 2 * (eg11 * eg22 + eg12 * eg21);
    }

    double compute(Ref<Vector> g, Ref<Matrix> h) const override
    {
        // actual configuration

        const HyperJet3D act_a1 = act_geometry<2>(1);
        const HyperJet3D act_a2 = act_geometry<2>(2);

        const HyperJet3D act_a1_1 = act_geometry<2>(3);
        const HyperJet3D act_a1_2 = act_geometry<2>(4);
        const HyperJet3D act_a2_2 = act_geometry<2>(5);

        HyperJet3D act_a3 = act_a1.cross(act_a2);
        const HyperJet act_da = act_a3.norm();
        act_a3 = act_a3 / act_da;

        const HyperJet act_a11 = act_a1.dot(act_a1);
        const HyperJet act_a22 = act_a2.dot(act_a2);
        const HyperJet act_a12 = act_a1.dot(act_a2);

        const HyperJet act_b11 = act_a1_1.dot(act_a3);
        const HyperJet act_b12 = act_a1_2.dot(act_a3);
        const HyperJet act_b22 = act_a2_2.dot(act_a3);

        const HyperJet3D delta_a(act_a11 - ref_a11, act_a22 - ref_a22, act_a12 - ref_a12);
        const HyperJet3D delta_b(ref_b11 - act_b11, ref_b12 - act_b12, ref_b22 - act_b22);

        const HyperJet3D eps = tm * delta_a * 0.5;
        const HyperJet3D kap = tm * delta_b;

        const HyperJet3D n = m_dm * eps;
        const HyperJet3D m = m_db * kap;

        const HyperJet p = (eps.dot(n) + kap.dot(m)) * 0.5 * m_weight;

        return hyperjet::explode(p, g, h);
    }

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Holder = Pointer<Type>;
        using Base = Objective;

        py::class_<Type, Base, Holder>(m, "IgaShell3PAD")
            .def(py::init<std::vector<Pointer<Point>>, Matrix, double, double, double, double>())
        ;
    }
};

} // namespace EQlib