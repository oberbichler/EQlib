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
        using namespace eqlib::iga_utilities;
        using Space = hyperjet::Space<2, double, -1>;
        using Scalar3hj = Space::Scalar;
        using Vector3hj = Space::Vector<3>;

        static_assert(0 <= TOrder && TOrder <= 2);

        auto result = Scalar3hj::zero(nb_variables());

        for (const auto& [shape_functions, ref_a, ref_b, transformation_matrix, weight] : m_data) {
            const Vector3hj act_a1 = evaluate_act_geometry_hj(m_nodes, shape_functions.row(1));
            const Vector3hj act_a2 = evaluate_act_geometry_hj(m_nodes, shape_functions.row(2));

            const Vector3hj act_a1_1 = evaluate_act_geometry_hj(m_nodes, shape_functions.row(3));
            const Vector3hj act_a1_2 = evaluate_act_geometry_hj(m_nodes, shape_functions.row(4));
            const Vector3hj act_a2_2 = evaluate_act_geometry_hj(m_nodes, shape_functions.row(5));

            const Vector3hj act_a3 = act_a1.cross(act_a2).normalized();

            const Vector3hj act_a(act_a1.dot(act_a1), act_a2.dot(act_a2), act_a1.dot(act_a2));
            const Vector3hj act_b(act_a1_1.dot(act_a3), act_a1_2.dot(act_a3), act_a2_2.dot(act_a3));

            const Vector3hj eps = transformation_matrix * (act_a - ref_a).transpose() / 2;
            const Vector3hj kap = transformation_matrix * (act_b - ref_b).transpose();

            result += (eps.dot(m_dm * eps.transpose()) + kap.dot(m_db * kap.transpose())) * weight;
        }

        g = result.g() / 2;
        h = result.h() / 2;
        return result.f() / 2;
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
}; // class Point

} // namespace eqlib