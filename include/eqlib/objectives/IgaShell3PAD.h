#pragma once

#include "IgaUtilities.h"

#include "../Node.h"
#include "../Objective.h"

#include <hyperjet/hyperjet.h>

namespace eqlib {

class IgaShell3PAD : public Objective {
private: // types
    using Type = IgaShell3PAD;

    struct Data {
        Matrix shape_functions;
        Eigen::Matrix<double, 1, 3> ref_a;
        Eigen::Matrix<double, 1, 3> ref_b;
        Eigen::Matrix3d transformation_matrix;
        double weight;
    };

private: // variables
    std::vector<Pointer<Node>> m_nodes;
    Eigen::Matrix3d m_dm;
    Eigen::Matrix3d m_db;
    std::vector<Data> m_data;

public: // constructor
    IgaShell3PAD(
        const std::vector<Pointer<Node>>& nodes,
        const double thickness,
        const double youngs_modulus,
        const double poissons_ratio)
        : m_nodes(nodes)
    {
        m_variables.reserve(length(nodes) * 3);

        for (const auto node : nodes) {
            m_variables.push_back(node->x());
            m_variables.push_back(node->y());
            m_variables.push_back(node->z());
        }

        m_dm << 1, poissons_ratio, 0, poissons_ratio, 1, 0, 0, 0, (1 - poissons_ratio) / 2;
        m_dm *= youngs_modulus * thickness / (1 - std::pow(poissons_ratio, 2));

        m_db << 1, poissons_ratio, 0, poissons_ratio, 1, 0, 0, 0, (1 - poissons_ratio) / 2;
        m_db *= youngs_modulus * std::pow(thickness, 3) / (12 * (1 - std::pow(poissons_ratio, 2)));
    }

public: // methods
    index add(const Matrix shape_functions, const double weight)
    {
        using namespace Eigen;
        using namespace eqlib::iga_utilities;

        const Vector3d ref_a1 = evaluate_ref_geometry(m_nodes, shape_functions.row(1));
        const Vector3d ref_a2 = evaluate_ref_geometry(m_nodes, shape_functions.row(2));

        const Vector3d ref_a1_1 = evaluate_ref_geometry(m_nodes, shape_functions.row(3));
        const Vector3d ref_a1_2 = evaluate_ref_geometry(m_nodes, shape_functions.row(4));
        const Vector3d ref_a2_2 = evaluate_ref_geometry(m_nodes, shape_functions.row(5));

        const Vector3d ref_a3 = ref_a1.cross(ref_a2).normalized();

        const Vector3d ref_a(ref_a1.dot(ref_a1), ref_a2.dot(ref_a2), ref_a1.dot(ref_a2));
        const Vector3d ref_b(ref_a1_1.dot(ref_a3), ref_a1_2.dot(ref_a3), ref_a2_2.dot(ref_a3));

        const Vector3d e1 = ref_a1.normalized();
        const Vector3d e2 = (ref_a2 - ref_a2.dot(e1) * e1).normalized();

        const double det = ref_a(0) * ref_a(1) - ref_a(2) * ref_a(2);

        const Vector3d g_ab_con(ref_a(1) / det, ref_a(0) / det, -ref_a(2) / det);

        const Vector3d g_con1 = g_ab_con[0] * ref_a1 + g_ab_con[2] * ref_a2;
        const Vector3d g_con2 = g_ab_con[2] * ref_a1 + g_ab_con[1] * ref_a2;

        const double eg11 = e1.dot(g_con1);
        const double eg12 = e1.dot(g_con2);
        const double eg21 = e2.dot(g_con1);
        const double eg22 = e2.dot(g_con2);

        Matrix3d transformation_matrix;
        transformation_matrix << eg11 * eg11, eg12 * eg12, 2 * eg11 * eg12, eg21 * eg21, eg22 * eg22, 2 * eg21 * eg22, 2 * eg11 * eg21, 2 * eg12 * eg22, 2 * (eg11 * eg22 + eg12 * eg21);

        m_data.emplace_back<Data>({shape_functions, ref_a, ref_b, transformation_matrix, weight});

        return m_data.size() - 1;
    }

    template <int TOrder>
    double compute(Ref<Vector> g, Ref<Matrix> h) const
    {
        using Space = hyperjet::Space<TOrder, double, 15>;

        static_assert(0 <= TOrder && TOrder <= 2);

        const index nb_nodes = length(m_nodes);

        Eigen::Matrix<double, Eigen::Dynamic, 3> locations(nb_nodes, 3);

        for (index i = 0; i < nb_nodes; i++) {
            locations.row(i) = m_nodes[i]->act_location();
        }

        double f = 0;

        if constexpr (TOrder > 0) {
            g.setZero();
        }

        if constexpr (TOrder > 1) {
            h.setZero();
        }

        for (const auto& [shape_functions, ref_a, ref_b, transformation_matrix, weight] : m_data) {
            const auto act_a1 = Space::template variables<0, 3>(shape_functions.row(1) * locations);
            const auto act_a2 = Space::template variables<3, 3>(shape_functions.row(2) * locations);

            const auto act_a1_1 = Space::template variables<6, 3>(shape_functions.row(3) * locations);
            const auto act_a1_2 = Space::template variables<9, 3>(shape_functions.row(4) * locations);
            const auto act_a2_2 = Space::template variables<12, 3>(shape_functions.row(5) * locations);

            const auto act_a3 = act_a1.cross(act_a2).normalized();

            const typename Space::template Vector<3> act_a(act_a1.dot(act_a1), act_a2.dot(act_a2), act_a1.dot(act_a2));
            const typename Space::template Vector<3> act_b(act_a1_1.dot(act_a3), act_a1_2.dot(act_a3), act_a2_2.dot(act_a3));

            const auto eps = transformation_matrix * (act_a - ref_a).transpose() * 0.5;
            const auto kap = transformation_matrix * (act_b - ref_b).transpose();

            const auto result = 0.5 * weight * (eps.dot(m_dm * eps) + kap.dot(m_db * kap));

            const index a = nb_nodes * 3;

            f += Space::f(result);

            if constexpr (TOrder > 0) {
                for (index r = 0; r < a; r++) {
                    const index rd = r % 3;
                    const index ri = r / 3;

                    for (index k = 0; k < 5; k++) {
                        g(r) += result.g(k * 3 + rd) * shape_functions(1 + k, ri);
                    }
                }
            }

            if constexpr (TOrder > 1) {
                for (index r = 0; r < a; r++) {
                    const index ri = r / 3;
                    const index rd = r % 3;

                    for (index s = r; s < a; s++) {
                        const index si = s / 3;
                        const index sd = s % 3;

                        for (index k = 0; k < 5; k++) {
                            for (index l = 0; l < 5; l++) {
                                h(r, s) += result.h(k * 3 + rd, l * 3 + sd) * shape_functions(1 + k, ri) * shape_functions(1 + l, si);
                            }
                        }
                    }
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

public: // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Holder = Pointer<Type>;
        using Base = Objective;

        py::class_<Type, Base, Holder>(m, "IgaShell3PAD")
            .def(py::init<std::vector<Pointer<Node>>, double, double, double>(), "nodes"_a, "thickness"_a, "youngs_modulus"_a, "poissons_ratio"_a)
            .def("add", &Type::add, "shape_functions"_a, "weight"_a);
    }
}; // class IgaShell3PAD

} // namespace eqlib