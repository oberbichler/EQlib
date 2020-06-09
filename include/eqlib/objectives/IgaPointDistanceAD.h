#pragma once

#include "IgaUtilities.h"

#include "../Node.h"
#include "../Objective.h"

#include <hyperjet/hyperjet.h>

namespace eqlib {

class IgaPointDistanceAD : public Objective {
private: // types
    using Type = IgaPointDistanceAD;

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
    IgaPointDistanceAD(
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
        using Space = hyperjet::Space<TOrder, double, 6>;

        static_assert(0 <= TOrder && TOrder <= 2);

        double f = 0;

        if constexpr (TOrder > 0) {
            g.setZero();
        }

        if constexpr (TOrder > 1) {
            h.setZero();
        }

        for (const auto& [shape_functions_a, shape_functions_b, weight] : m_data) {
            const auto x_a = Space::variables<0, 3>(evaluate_act_geometry(m_nodes_a, shape_functions_a.row(0)));
            const auto x_b = Space::variables<3, 3>(evaluate_act_geometry(m_nodes_b, shape_functions_b.row(0)));

            const auto delta = x_a - x_b;

            const auto result = 0.5 * weight * delta.squaredNorm();

            f += Space::f(result);

            const index a = shape_functions_a.cols() * 3;
            const index b = shape_functions_b.cols() * 3;

            if constexpr (TOrder > 0) {
                for (index r = 0; r < a; r++) {
                    const index rd = r % 3;
                    const index ri = r / 3;

                    g(0 + r) = result.g(0 + rd) * shape_functions_a(0, ri);
                }

                for (index r = 0; r < b; r++) {
                    const index rd = r % 3;
                    const index ri = r / 3;

                    g(a + r) = result.g(3 + rd) * shape_functions_b(0, ri);
                }
            }

            if constexpr (TOrder > 1) {
                for (index r = 0; r < a; r++) {
                    const index rd = r % 3;
                    const index ri = r / 3;

                    for (index s = r; s < a; s++) {
                        const index sd = s % 3;
                        const index si = s / 3;

                        if (rd != sd) {
                            continue;
                        }

                        h(0 + r, 0 + s) = result.h(0 + rd, 0 + sd) * shape_functions_a(0, ri) * shape_functions_a(0, si);
                    }

                    for (index s = 0; s < b; s++) {
                        const index sd = s % 3;
                        const index si = s / 3;

                        if (rd != sd) {
                            continue;
                        }

                        h(0 + r, a + s) = result.h(0 + rd, 3 + sd) * shape_functions_a(0, ri) * shape_functions_b(0, si);
                    }
                }

                for (index r = 0; r < b; r++) {
                    const index rd = r % 3;
                    const index ri = r / 3;

                    for (index s = r; s < b; s++) {
                        const index sd = s % 3;
                        const index si = s / 3;

                        if (rd != sd) {
                            continue;
                        }

                        h(a + r, a + s) = result.h(3 + rd, 3 + sd) * shape_functions_b(0, ri) * shape_functions_b(0, si);
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

        py::class_<Type, Base, Holder>(m, "IgaPointDistanceAD")
            .def(py::init<std::vector<Pointer<Node>>, std::vector<Pointer<Node>>>(), "nodes_a"_a, "nodes_b"_a)
            .def("add", &Type::add, "shape_functions"_a, "shape_functions_b"_a, "weight"_a);
    }
}; // class Point

} // namespace eqlib