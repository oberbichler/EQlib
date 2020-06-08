#pragma once

#include "IgaUtilities.h"

#include "../Node.h"
#include "../Objective.h"

#include <hyperjet/hyperjet.h>

namespace eqlib {

class IgaNormalDistanceAD : public Objective {
private: // types
    using Type = IgaNormalDistanceAD;

    struct Data {
        Matrix shape_functions_a;
        Matrix shape_functions_b;
        double weight;
    };

private: // variables
    std::vector<Pointer<Node>> m_nodes_a;
    std::vector<Pointer<Node>> m_nodes_b;
    std::vector<Data> m_data;

public: // constructor
    IgaNormalDistanceAD(
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
    index add(const Matrix shape_functions_a, const Matrix shape_functions_b, const double weight)
    {
        m_data.emplace_back<Data>({shape_functions_a, shape_functions_b, weight});
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

        for (const auto& [shape_functions_a, shape_functions_b, weight] : m_data) {
            const auto a1_a = evaluate_act_geometry_hj_a(m_nodes_a, shape_functions_a.row(1), nb_variables());
            const auto a2_a = evaluate_act_geometry_hj_a(m_nodes_a, shape_functions_a.row(2), nb_variables());

            const auto a3_a = a1_a.cross(a2_a).normalized();

            const auto a1_b = evaluate_act_geometry_hj_b(m_nodes_b, shape_functions_b.row(1), nb_variables());
            const auto a2_b = evaluate_act_geometry_hj_b(m_nodes_b, shape_functions_b.row(2), nb_variables());

            const auto a3_b = a1_b.cross(a2_b).normalized();

            const auto delta = a3_b - a3_a;

            result += delta.dot(delta) * weight;
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

        py::class_<Type, Base, Holder>(m, "IgaNormalDistanceAD")
            .def(py::init<std::vector<Pointer<Node>>, std::vector<Pointer<Node>>>(), "nodes_a"_a, "nodes_b"_a)
            .def("add", &Type::add, "shape_functions"_a, "shape_functions_b"_a, "weight"_a);
    }
}; // class Point

} // namespace eqlib