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
        using Space = hyperjet::Space<TOrder, double, 12>;

        static_assert(0 <= TOrder && TOrder <= 2);

        double f = 0;

        for (const auto& [shape_functions_a, shape_functions_b, weight] : m_data) {
            const auto a1_a = Space::template variables<0, 3>(evaluate_act_geometry(m_nodes_a, shape_functions_a.row(1)));
            const auto a2_a = Space::template variables<3, 3>(evaluate_act_geometry(m_nodes_a, shape_functions_a.row(2)));
            
            const auto a1_b = Space::template variables<6, 3>(evaluate_act_geometry(m_nodes_b, shape_functions_b.row(1)));
            const auto a2_b = Space::template variables<9, 3>(evaluate_act_geometry(m_nodes_b, shape_functions_b.row(2)));

            const auto a3_a = a1_a.cross(a2_a).normalized();
            const auto a3_b = a1_b.cross(a2_b).normalized();

            const auto delta = a3_b - a3_a;

            const auto result = 0.5 * weight * delta.dot(delta);
            
            const index a = shape_functions_a.cols() * 3;
            const index b = shape_functions_b.cols() * 3;

            if constexpr(TOrder > 0) {
                for (index r = 0; r < a; r++) {
                    const index rd = r % 3;
                    const index ri = r / 3;

                    g(0 + r) = result.g(0 + rd) * shape_functions_a(1, ri) + result.g(3 + rd) * shape_functions_a(2, ri);
                }

                for (index r = 0; r < b; r++) {
                    const index rd = r % 3;
                    const index ri = r / 3;

                    g(a + r) = result.g(6 + rd) * shape_functions_b(1, ri) + result.g(9 + rd) * shape_functions_b(2, ri);
                }
            }

            if constexpr(TOrder > 1) {
                for (index r = 0; r < a; r++) {
                    const index rd = r % 3;
                    const index ri = r / 3;

                    for (index s = r; s < a; s++) {
                        const index sd = s % 3;
                        const index si = s / 3;

                        h(0 + r, 0 + s) = result.h(0 + rd, 0 + sd) * shape_functions_a(1, ri) * shape_functions_a(1, si)
                                        + result.h(3 + rd, 0 + sd) * shape_functions_a(2, ri) * shape_functions_a(1, si)
                                        + result.h(0 + rd, 3 + sd) * shape_functions_a(1, ri) * shape_functions_a(2, si)
                                        + result.h(3 + rd, 3 + sd) * shape_functions_a(2, ri) * shape_functions_a(2, si);
                    }

                    for (index s = 0; s < b; s++) {
                        const index sd = s % 3;
                        const index si = s / 3;

                        h(0 + r, a + s) = result.h(0 + rd, 6 + sd) * shape_functions_a(1, ri) * shape_functions_b(1, si)
                                        + result.h(3 + rd, 6 + sd) * shape_functions_a(2, ri) * shape_functions_b(1, si)
                                        + result.h(0 + rd, 9 + sd) * shape_functions_a(1, ri) * shape_functions_b(2, si)
                                        + result.h(3 + rd, 9 + sd) * shape_functions_a(2, ri) * shape_functions_b(2, si);
                    }
                }

                for (index r = 0; r < b; r++) {
                    const index rd = r % 3;
                    const index ri = r / 3;

                    for (index s = r; s < b; s++) {
                        const index sd = s % 3;
                        const index si = s / 3;

                        h(a + r, a + s) = result.h(6 + rd, 6 + sd) * shape_functions_b(1, ri) * shape_functions_b(1, si)
                                        + result.h(9 + rd, 6 + sd) * shape_functions_b(2, ri) * shape_functions_b(1, si)
                                        + result.h(6 + rd, 9 + sd) * shape_functions_b(1, ri) * shape_functions_b(2, si)
                                        + result.h(9 + rd, 9 + sd) * shape_functions_b(2, ri) * shape_functions_b(2, si);
                    }
                }
            }

            f += Space::f(result);
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

        py::class_<Type, Base, Holder>(m, "IgaNormalDistanceAD")
            .def(py::init<std::vector<Pointer<Node>>, std::vector<Pointer<Node>>>(), "nodes_a"_a, "nodes_b"_a)
            .def("add", &Type::add, "shape_functions"_a, "shape_functions_b"_a, "weight"_a);
    }
}; // class IgaNormalDistanceAD

} // namespace eqlib