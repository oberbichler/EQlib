#pragma once

#include "IgaUtilities.h"

#include "../Node.h"
#include "../Objective.h"

namespace eqlib {

class IgaPointLocation : public Objective {
private: // types
    using Type = IgaPointLocation;

    struct Data {
        Matrix shape_functions;
        Vector3D location;
        double weight;
    };

private: // variables
    std::vector<Pointer<Node>> m_nodes;
    std::vector<Data> m_data;

public: // constructor
    IgaPointLocation(std::vector<Pointer<Node>> nodes)
        : m_nodes(nodes)
    {
        m_variables.resize(length(nodes) * 3);
        for (index i = 0; i < length(nodes); i++) {
            m_variables[i * 3 + 0] = nodes[i]->x();
            m_variables[i * 3 + 1] = nodes[i]->y();
            m_variables[i * 3 + 2] = nodes[i]->z();
        }
    }

public: // methods
    index add(const Matrix shape_functions, const Vector target, const double weight)
    {
        m_data.emplace_back<Data>({shape_functions, target, weight});
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

        for (const auto& [shape_functions, target, weight] : m_data) {
            const Vector3D act_x = evaluate_act_geometry(m_nodes, shape_functions.row(0));

            const Vector3D delta = act_x - target;

            f += 0.5 * weight * delta.dot(delta);

            if constexpr (TOrder > 0) {
                for (index i = 0; i < length(m_nodes); i++) {
                    g(i * 3 + 0) += delta[0] * shape_functions(0, i) * weight;
                    g(i * 3 + 1) += delta[1] * shape_functions(0, i) * weight;
                    g(i * 3 + 2) += delta[2] * shape_functions(0, i) * weight;
                }
            }

            if constexpr (TOrder > 1) {
                for (index i = 0; i < length(m_nodes); i++) {
                    for (index j = 0; j < length(m_nodes); j++) {
                        h(i * 3 + 0, j * 3 + 0) += shape_functions(0, i) * shape_functions(0, j) * weight;
                        h(i * 3 + 1, j * 3 + 1) += shape_functions(0, i) * shape_functions(0, j) * weight;
                        h(i * 3 + 2, j * 3 + 2) += shape_functions(0, i) * shape_functions(0, j) * weight;
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

        py::class_<Type, Base, Holder>(m, "IgaPointLocation")
            .def(py::init<std::vector<Pointer<Node>>>(), "nodes"_a)
            .def("add", &Type::add, "shape_functions"_a, "target"_a, "weight"_a);
    }
}; // class IgaPointLocation

} // namespace eqlib