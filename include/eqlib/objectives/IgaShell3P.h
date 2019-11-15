#pragma once

#include <Eigen/Geometry>

#include <eqlib/Node.h>
#include <eqlib/Objective.h>
#include <eqlib/Variable.h>

#include <vector>

namespace eqlib {

class IgaShell3P : public Objective
{
private:    // types
    using Type = IgaShell3P;
    using InputData = std::tuple<Matrix, double, double, double, double>;

    struct Configuration
    {
        Configuration(const Vector3D& a1, const Vector3D& a2, const Vector3D& a1_1, const Vector3D& a1_2, const Vector3D& a2_2)
        {
            a11 = a1.dot(a1);
            a12 = a1.dot(a2);
            a22 = a2.dot(a2);

            a1_x_a2 = a1.cross(a2);
            da = a1.cross(a2).norm();

            a3 = a1_x_a2 / da;

            b22 = a1_1.dot(a3);
            b22 = a1_2.dot(a3);
            b22 = a2_2.dot(a3);

            const Vector3D e1 = a1.normalized();
            const Vector3D e2 = (a2 - a2.dot(e1) * e1).normalized();

            const double det = a11 * a22 - a12 * a12;

            const Vector3D g_ab_con(a22 / det, a11 / det, -a12 / det);

            const Vector3D g_con1 = g_ab_con[0] * a1 + g_ab_con[2] * a2;
            const Vector3D g_con2 = g_ab_con[2] * a1 + g_ab_con[1] * a2;

            const double eg11 = e1.dot(g_con1);
            const double eg12 = e1.dot(g_con2);
            const double eg21 = e2.dot(g_con1);
            const double eg22 = e2.dot(g_con2);

            tm(0, 0) = eg11 * eg11;
            tm(0, 1) = eg12 * eg12;
            tm(0, 2) = 2 * eg11 * eg12;
            tm(1, 0) = eg21 * eg21;
            tm(1, 1) = eg22 * eg22;
            tm(1, 2) = 2 * eg21 * eg22;
            tm(2, 0) = 2 * eg11 * eg21;
            tm(2, 1) = 2 * eg12 * eg22;
            tm(2, 2) = 2 * (eg11 * eg22 + eg12 * eg21);
        }

        double a11;
        double a12;
        double a22;

        double b11;
        double b12;
        double b22;

        Vector3D a1_x_a2;
        double da;
        Vector3D a3;

        Eigen::Matrix3d tm;
    };

    struct Data
    {
        Matrix shape_functions;
        Configuration ref;
        Eigen::Matrix3d dm;
        Eigen::Matrix3d db;
        double weight;
    };

private:    // variables
    std::vector<Pointer<Node>> m_nodes;
    std::vector<Data> m_data;
    double m_weight;

    Vector3D ref_geometry(index i, Ref<const Matrix> shape_functions) const
    {
        Vector3D value = Vector3D::Zero();

        for (index j = 0; j < length(m_nodes); j++) {
            value += m_nodes[j]->ref_location() * shape_functions(i, j);
        }

        return value;
    }

    Vector3D act_geometry(index i, Ref<const Matrix> shape_functions) const
    {
        Vector3D value = Vector3D::Zero();

        for (index j = 0; j < length(m_nodes); j++) {
            value += m_nodes[j]->act_location() * shape_functions(i, j);
        }

        return value;
    }

public:     // constructor
    IgaShell3P(
        std::vector<Pointer<Node>> nodes,
        std::vector<InputData> data)
    : m_nodes(nodes)
    {
        m_variables.resize(length(nodes) * 3);
        for (index i = 0; i < length(nodes); i++) {
            m_variables[i * 3 + 0] = nodes[i]->x();
            m_variables[i * 3 + 1] = nodes[i]->y();
            m_variables[i * 3 + 2] = nodes[i]->z();
        }

        for (const auto& [shape_functions, thickness, youngs_modulus, poisson_ratio, weight] : data) {
            const Vector3D ref_a1 = ref_geometry(1, shape_functions);
            const Vector3D ref_a2 = ref_geometry(2, shape_functions);
            const Vector3D ref_a1_1 = ref_geometry(3, shape_functions);
            const Vector3D ref_a1_2 = ref_geometry(4, shape_functions);
            const Vector3D ref_a2_2 = ref_geometry(5, shape_functions);

            const Vector3D ref_a3 = ref_a1.cross(ref_a2).normalized();

            Configuration ref{ref_a1, ref_a2, ref_a1_1, ref_a1_2, ref_a2_2};

            Eigen::Matrix3d dm;
            dm << 1.0, poisson_ratio, 0,
                  poisson_ratio, 1.0, 0,
                  0, 0, (1.0 - poisson_ratio) / 2.0;
            dm *= youngs_modulus * thickness / (1.0 - poisson_ratio * poisson_ratio);

            Eigen::Matrix3d db;
            db << 1.0, poisson_ratio, 0,
                  poisson_ratio, 1.0, 0,
                  0, 0, (1.0 - poisson_ratio) / 2.0;
            db *= youngs_modulus * std::pow(thickness, 3) / (12.0 * (1.0 - poisson_ratio * poisson_ratio));

            m_data.push_back({shape_functions, ref, dm, db, weight});
        }
    }

    template <int TOrder>
    double compute(Ref<Vector> g, Ref<Matrix> h) const
    {
        static_assert(0 <= TOrder && TOrder <= 2);

        double f = 0;
        g.setZero();
        h.setZero();

        const index nb_nodes = length(m_nodes);
        const index nb_dofs = nb_nodes * 3;

        Vector s_g3dg3(nb_dofs);
        Vector s_g3dg3lg3_3(nb_dofs);

        std::vector<Vector3D> s_dg3(nb_dofs);
        std::vector<Vector3D> s_de_ca(nb_dofs);
        std::vector<Vector3D> s_dk_ca(nb_dofs);
        std::vector<Vector3D> s_dn(nb_dofs);

        std::vector<Vector3D> s_dde_ca(nb_dofs * nb_dofs);
        std::vector<Vector3D> s_ddk_ca(nb_dofs * nb_dofs);

        for (const auto& [shape_functions, ref, dm, db, weight] : m_data) {
            const Vector3D act_a1 = act_geometry(1, shape_functions);
            const Vector3D act_a2 = act_geometry(2, shape_functions);
            const Vector3D act_a1_1 = act_geometry(3, shape_functions);
            const Vector3D act_a1_2 = act_geometry(4, shape_functions);
            const Vector3D act_a2_2 = act_geometry(5, shape_functions);

            const Configuration act(act_a1, act_a2, act_a1_1, act_a1_2, act_a2_2);

            const Vector3D e_cu((act.a11 - ref.a11) * 0.5, (act.a22 - ref.a22) * 0.5, (act.a12 - ref.a12) * 0.5);
            const Vector3D k_cu(act.b11 - ref.b11, act.b22 - ref.b22, act.b12 - ref.b12);

            const Vector3D e_ca = ref.tm * e_cu;
            const Vector3D k_ca = ref.tm * k_cu;
    
            for (index r = 0; r < nb_dofs; r++) {
                const index shape_function_index_r = r / 3;
                const index dof_type_index_r = r % 3;

                Vector3D S_dg_1 = Vector3D::Zero();
                Vector3D S_dg_2 = Vector3D::Zero();

                S_dg_1[dof_type_index_r] = shape_functions(1, shape_function_index_r);
                S_dg_2[dof_type_index_r] = shape_functions(2, shape_function_index_r);

                // strain
                Vector3D dE_cu;
                dE_cu[0] = shape_functions(1, shape_function_index_r) * act_a1[dof_type_index_r];
                dE_cu[1] = shape_functions(2, shape_function_index_r) * act_a2[dof_type_index_r];
                dE_cu[2] = 0.5 * (shape_functions(1, shape_function_index_r) * act_a2[dof_type_index_r] + shape_functions(2, shape_function_index_r) * act_a1[dof_type_index_r]);

                s_de_ca[r] = ref.tm * dE_cu;

                // curvature
                s_dg3[r] = S_dg_1.cross(act_a2) + act_a1.cross(S_dg_2);

                s_g3dg3[r] = ref.a1_x_a2.dot(s_dg3[r]);
                s_g3dg3lg3_3[r] = s_g3dg3[r] / std::pow(ref.da, 3);

                s_dn[r] = s_dg3[r] / ref.da - ref.a1_x_a2 * s_g3dg3lg3_3[r];

                Vector3D dK_cu;
                dK_cu[0] = shape_functions(3, shape_function_index_r) * act.a3[dof_type_index_r] + act_a1_1.dot(s_dn[r]);
                dK_cu[1] = shape_functions(4, shape_function_index_r) * act.a3[dof_type_index_r] + act_a1_2.dot(s_dn[r]);
                dK_cu[2] = shape_functions(5, shape_function_index_r) * act.a3[dof_type_index_r] + act_a2_2.dot(s_dn[r]);

                s_dk_ca[r] = ref.tm * dK_cu;
            }

            for (index r = 0; r < nb_dofs; r++) {
                const index shape_function_index_r = r / 3;
                const index dof_type_index_r = r % 3;

                for (index s = 0; s <= r; s++) {
                    const index shape_index_s = s / 3;
                    const index dof_type_index_s = s % 3;

                    // strain

                    if (dof_type_index_r == dof_type_index_s) {
                        Vector3D ddE_cu;

                        ddE_cu[0] = shape_functions(1, shape_function_index_r) * shape_functions(1, shape_index_s);
                        ddE_cu[1] = shape_functions(2, shape_function_index_r) * shape_functions(2, shape_index_s);
                        ddE_cu[2] = 0.5 * (shape_functions(1, shape_function_index_r) * shape_functions(2, shape_index_s) + shape_functions(2, shape_function_index_r) * shape_functions(1, shape_index_s));

                        s_dde_ca[r * nb_dofs + s] = ref.tm * ddE_cu;
                    } else {
                        s_dde_ca[r * nb_dofs + s].setZero();
                    }

                    // curvature

                    Vector3D ddg3 = Vector3D::Zero();

                    if (dof_type_index_r != dof_type_index_s) {
                        const index ddg3_i = 3 - dof_type_index_r - dof_type_index_s;

                        const double ddg3_value = shape_functions(1, shape_function_index_r) * shape_functions(2, shape_index_s) - shape_functions(1, shape_index_s) * shape_functions(2, shape_function_index_r);

                        if ((dof_type_index_s == dof_type_index_r + 1) || (dof_type_index_r == dof_type_index_s + 2)) {
                            ddg3[ddg3_i] = ddg3_value;
                        } else {
                            ddg3[ddg3_i] = -ddg3_value;
                        }
                    }

                    const double c = -(ddg3.dot(ref.a1_x_a2) + s_dg3[r].dot(s_dg3[s])) / std::pow(ref.da, 3);

                    const double d = 3.0 * s_g3dg3[r] * s_g3dg3[s] / std::pow(ref.da, 5);

                    const Vector3D ddn = ddg3 / ref.da - s_g3dg3lg3_3[s] * s_dg3[r] - s_g3dg3lg3_3[r] * s_dg3[s] + (c + d) * ref.a1_x_a2;

                    Vector3D ddK_cu;
                    ddK_cu[0] = shape_functions(3, shape_function_index_r) * s_dn[s][dof_type_index_r] + shape_functions(3, shape_index_s) * s_dn[r][dof_type_index_s] + act_a1_1[0] * ddn[0] + act_a1_1[1] * ddn[1] + act_a1_1[2] * ddn[2];
                    ddK_cu[1] = shape_functions(4, shape_function_index_r) * s_dn[s][dof_type_index_r] + shape_functions(4, shape_index_s) * s_dn[r][dof_type_index_s] + act_a1_2[0] * ddn[0] + act_a1_2[1] * ddn[1] + act_a1_2[2] * ddn[2];
                    ddK_cu[2] = shape_functions(5, shape_function_index_r) * s_dn[s][dof_type_index_r] + shape_functions(5, shape_index_s) * s_dn[r][dof_type_index_s] + act_a2_2[0] * ddn[0] + act_a2_2[1] * ddn[1] + act_a2_2[2] * ddn[2];

                    s_ddk_ca[r * nb_dofs + s] = ref.tm * ddK_cu;
                }
            }

            const Vector3D N_ca = dm * e_ca;
            const Vector3D M_ca = db * k_ca;

            for(index r = 0; r < nb_dofs; r++) {
                const Vector3D dN_ca = dm * s_de_ca[r];
                const Vector3D dM_ca = db * s_dk_ca[r];

                if constexpr(TOrder > 1) {
                    for (index s = 0; s < nb_dofs; s++) {
                        // membrane stiffness
                        const double S_kem = dN_ca.dot(s_de_ca[s]) + N_ca.dot(s_dde_ca[r * nb_dofs + s]);

                        // bending stiffness
                        const double S_keb = dM_ca.dot(s_dk_ca[s]) + M_ca.dot(s_ddk_ca[r * nb_dofs + s]);

                        h(r, s) = (S_kem + S_keb) * weight;

                        // symmetry
                        h(s, r) = h(r, s);
                    }
                }

                if constexpr(TOrder > 0) {
                    g[r] = weight * (N_ca.dot(s_de_ca[r]) + M_ca.dot(s_dk_ca[r]));
                }
            }
        }

        return f;
    }

    double compute(Ref<Vector> g, Ref<Matrix> h) const override
    {
        if (g.size() == 0) {
            return compute<0>(g, h);
        } else if (h.size() == 0) {
            return compute<1>(g, h);
        } else {
            return compute<2>(g, h);
        }
    }

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Holder = Pointer<Type>;
        using Base = Objective;

        py::class_<Type, Base, Holder>(m, "IgaShell3P")
            .def(py::init<std::vector<Pointer<Node>>, std::vector<InputData>>())
        ;
    }
};

} // namespace eqlib