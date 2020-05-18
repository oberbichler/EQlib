#pragma once

#include "IgaUtilities.h"

#include "../Node.h"
#include "../Objective.h"

#include <hyperjet/hyperjet.h>

namespace eqlib {

class IgaPointDistance : public Objective {
private: // types
    using Type = IgaPointDistance;

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
    IgaPointDistance(
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

        static_assert(0 <= TOrder && TOrder <= 2);

        double f = 0;
        g.setZero();
        h.setZero();

        for (const auto& [shape_functions_a, shape_functions_b, weight] : m_data) {
            Vector3D point_a = evaluate_act_geometry(m_nodes_a, shape_functions_a.row(0));
            Vector3D point_b = evaluate_act_geometry(m_nodes_b, shape_functions_b.row(0));

            const Vector3D delta = point_b - point_a;

            f += delta.squaredNorm() * weight / 2;

            if constexpr (TOrder > 0) {
                index value_index = 0;

                for (index i = 0; i < length(m_nodes_a); i++) {
                    g(value_index++) -= delta[0] * shape_functions_a(0, i) * weight;
                    g(value_index++) -= delta[1] * shape_functions_a(0, i) * weight;
                    g(value_index++) -= delta[2] * shape_functions_a(0, i) * weight;
                }

                for (index i = 0; i < length(m_nodes_b); i++) {
                    g(value_index++) += delta[0] * shape_functions_b(0, i) * weight;
                    g(value_index++) += delta[1] * shape_functions_b(0, i) * weight;
                    g(value_index++) += delta[2] * shape_functions_b(0, i) * weight;
                }
            }

            if constexpr (TOrder > 1) {
                index row_index = 0;

                for (index i = 0; i < length(m_nodes_a); i++) {
                    index col_index = 0;

                    for (index j = 0; j < length(m_nodes_a); j++) {
                        h(row_index + 0, col_index++) += shape_functions_a(0, i) * shape_functions_a(0, j) * weight;
                        h(row_index + 1, col_index++) += shape_functions_a(0, i) * shape_functions_a(0, j) * weight;
                        h(row_index + 2, col_index++) += shape_functions_a(0, i) * shape_functions_a(0, j) * weight;
                    }

                    for (index j = 0; j < length(m_nodes_b); j++) {
                        h(row_index + 0, col_index++) -= shape_functions_a(0, i) * shape_functions_b(0, j) * weight;
                        h(row_index + 1, col_index++) -= shape_functions_a(0, i) * shape_functions_b(0, j) * weight;
                        h(row_index + 2, col_index++) -= shape_functions_a(0, i) * shape_functions_b(0, j) * weight;
                    }

                    row_index += 3;
                }

                for (index i = 0; i < length(m_nodes_b); i++) {
                    index col_index = 0;

                    for (index j = 0; j < length(m_nodes_a); j++) {
                        h(row_index + 0, col_index++) -= shape_functions_b(0, i) * shape_functions_a(0, j) * weight;
                        h(row_index + 1, col_index++) -= shape_functions_b(0, i) * shape_functions_a(0, j) * weight;
                        h(row_index + 2, col_index++) -= shape_functions_b(0, i) * shape_functions_a(0, j) * weight;
                    }

                    for (index j = 0; j < length(m_nodes_b); j++) {
                        h(row_index + 0, col_index++) += shape_functions_b(0, i) * shape_functions_b(0, j) * weight;
                        h(row_index + 1, col_index++) += shape_functions_b(0, i) * shape_functions_b(0, j) * weight;
                        h(row_index + 2, col_index++) += shape_functions_b(0, i) * shape_functions_b(0, j) * weight;
                    }

                    row_index += 3;
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

        py::class_<Type, Base, Holder>(m, "IgaPointDistance")
            .def(py::init<std::vector<Pointer<Node>>, std::vector<Pointer<Node>>>(), "nodes_a"_a, "nodes_b"_a)
            .def("add", &Type::add, "shape_functions"_a, "shape_functions_b"_a, "weight"_a);
    }
}; // class Point

} // namespace eqlib