#pragma once

#include "IgaUtilities.h"

#include "../Node.h"
#include "../Objective.h"

#include <hyperjet/hyperjet.h>

namespace eqlib {

class IgaRotationCouplingAD : public Objective {
private: // types
    using Type = IgaRotationCouplingAD;

    struct Data {
        Matrix shape_functions_a;
        Matrix shape_functions_b;
        Vector3D ref_a3_a;
        Vector3D ref_a3_b;
        Vector3D axis;
        double weight;
    };

private: // variables
    std::vector<Pointer<Node>> m_nodes_a;
    std::vector<Pointer<Node>> m_nodes_b;
    std::vector<Data> m_data;

public: // constructor
    IgaRotationCouplingAD(
        std::vector<Pointer<Node>> nodes_a,
        std::vector<Pointer<Node>> nodes_b)
        : m_nodes_a(nodes_a)
        , m_nodes_b(nodes_b)
    {
        m_variables.reserve(length(nodes_a) * 3 + length(nodes_b) * 3);

        for (const auto node : nodes_a) {
            m_variables.push_back(node->x());
            m_variables.push_back(node->y());
            m_variables.push_back(node->z());
        }

        for (const auto node : nodes_b) {
            m_variables.push_back(node->x());
            m_variables.push_back(node->y());
            m_variables.push_back(node->z());
        }
    }

public: // methods
    index add(const Matrix shape_functions_a, const Matrix shape_functions_b, const Vector3D axis, const double weight)
    {
        using namespace iga_utilities;

        const Vector3D ref_a1_a = evaluate_ref_geometry(m_nodes_a, shape_functions_a.row(1));
        const Vector3D ref_a2_a = evaluate_ref_geometry(m_nodes_a, shape_functions_a.row(2));
        const Vector3D ref_a3_a = ref_a1_a.cross(ref_a2_a).normalized();

        const Vector3D ref_a1_b = evaluate_ref_geometry(m_nodes_b, shape_functions_b.row(1));
        const Vector3D ref_a2_b = evaluate_ref_geometry(m_nodes_b, shape_functions_b.row(2));
        const Vector3D ref_a3_b = ref_a1_b.cross(ref_a2_b).normalized();

        const Vector3D unit_axis = axis.normalized();

        m_data.emplace_back<Data>({shape_functions_a, shape_functions_b, ref_a3_a, ref_a3_b, unit_axis, weight});
        return m_data.size() - 1;
    }

    template <int TOrder>
    double compute(Ref<Vector> g, Ref<Matrix> h) const
    {
        using namespace eqlib::iga_utilities;
        using Space = hyperjet::Space<2, double, -1>;

        static_assert(0 <= TOrder && TOrder <= 2);

        auto result = Space::Scalar::zero(nb_variables());

        for (const auto& [shape_functions_a, shape_functions_b, ref_a3_a, ref_a3_b, axis, weight] : m_data) {
            const auto a1_a = evaluate_act_geometry_hj_a(m_nodes_a, shape_functions_a.row(1), nb_variables());
            const auto a2_a = evaluate_act_geometry_hj_a(m_nodes_a, shape_functions_a.row(2), nb_variables());

            const auto a3_a = a1_a.cross(a2_a).normalized();

            const auto a1_b = evaluate_act_geometry_hj_b(m_nodes_b, shape_functions_b.row(1), nb_variables());
            const auto a2_b = evaluate_act_geometry_hj_b(m_nodes_b, shape_functions_b.row(2), nb_variables());

            const auto a3_b = a1_b.cross(a2_b).normalized();

            const auto w_a = a3_a - ref_a3_a.transpose();
            const auto w_b = a3_b - ref_a3_b.transpose();

            const auto omega_a = ref_a3_a.cross(w_a);
            const auto omega_b = ref_a3_b.cross(w_b);

            const auto angle_a = omega_a.dot(axis).asin();
            const auto angle_b = omega_b.dot(axis).asin();

            const auto angular_difference = angle_a - angle_b;

            result += angular_difference.pow(2) * weight;
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

        py::class_<Type, Base, Holder>(m, "IgaRotationCouplingAD")
            .def(py::init<std::vector<Pointer<Node>>, std::vector<Pointer<Node>>>(), "nodes_a"_a, "nodes_b"_a)
            .def("add", &Type::add, "shape_functions"_a, "shape_functions_b"_a, "axis"_a, "weight"_a);
    }
}; // class IgaRotationCouplingAD

} // namespace eqlib