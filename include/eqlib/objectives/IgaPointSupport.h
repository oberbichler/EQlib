#pragma once

#include <eqlib/Node.h>
#include <eqlib/Objective.h>
#include <eqlib/Variable.h>

#include <vector>

namespace eqlib {

class IgaPointSupport : public Objective
{
private:    // types
    using Type = IgaPointSupport;
    using Data = std::tuple<Matrix, Vector3D, double>;

private:    // variables
    std::vector<Pointer<Node>> m_nodes;
    std::vector<Data> m_data;

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
    IgaPointSupport(
        std::vector<Pointer<Node>> nodes,
        std::vector<Data> data)
    : m_nodes(nodes)
    , m_data(data)
    {
        m_variables.resize(length(nodes) * 3);
        for (index i = 0; i < length(nodes); i++) {
            m_variables[i * 3 + 0] = nodes[i]->x();
            m_variables[i * 3 + 1] = nodes[i]->y();
            m_variables[i * 3 + 2] = nodes[i]->z();
        }
    }

    template <int TOrder>
    double compute(Ref<Vector> g, Ref<Matrix> h) const
    {
        static_assert(0 <= TOrder && TOrder <= 2);

        double f = 0;
        g.setZero();
        h.setZero();

        for (const auto& [shape_functions, target, weight] : m_data) {
            const Vector3D act_x = act_geometry(0, shape_functions);

            const Vector3D delta = act_x - target;

            f += delta.dot(delta) * weight / 2;

            if constexpr(TOrder > 0) {
                for (index i = 0; i < length(m_nodes); i++) {
                    g(i * 3 + 0) += delta[0] * shape_functions(0, i) * weight;
                    g(i * 3 + 1) += delta[1] * shape_functions(0, i) * weight;
                    g(i * 3 + 2) += delta[2] * shape_functions(0, i) * weight;
                }
            }

            if constexpr(TOrder > 1) {
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

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Holder = Pointer<Type>;
        using Base = Objective;

        py::class_<Type, Base, Holder>(m, "IgaPointSupport")
            .def(py::init<std::vector<Pointer<Node>>, std::vector<Data>>(), "nodes"_a, "data"_a)
        ;
    }
};

} // namespace eqlib