#pragma once

#include <Eigen/Geometry>

#include <eqlib/Objective.h>
#include <eqlib/Node.h>
#include <eqlib/Variable.h>

#include <hyperjet/hyperjet.h>

#include <vector>

namespace eqlib {

class IgaDisplacementCoupling : public Objective
{
private:    // types
    using Type = IgaDisplacementCoupling;

    using Data = std::tuple<Matrix, Matrix, double>;

private:    // variables
    std::vector<Pointer<Node>> m_nodes_a;
    std::vector<Pointer<Node>> m_nodes_b;
    std::vector<Data> m_data;

    Vector3D ref_geometry(const std::vector<Pointer<Node>>& nodes, Ref<const Matrix> shape_functions, const index i) const
    {
        Vector3D value = Vector3D::Zero();

        for (index j = 0; j < length(nodes); j++) {
            value += nodes[j]->ref_location() * shape_functions(i, j);
        }

        return value;
    }

    Vector3D act_geometry(const std::vector<Pointer<Node>>& nodes, Ref<const Matrix> shape_functions, const index i) const
    {
        Vector3D value = Vector3D::Zero();

        for (index j = 0; j < length(nodes); j++) {
            value += nodes[j]->act_location() * shape_functions(i, j);
        }

        return value;
    }

public:     // constructor
    IgaDisplacementCoupling(
        std::vector<Pointer<Node>> nodes_a,
        std::vector<Pointer<Node>> nodes_b,
        std::vector<Data> data)
    : m_nodes_a(nodes_a)
    , m_nodes_b(nodes_b)
    , m_data(data)
    {
        const index nb_dofs = (length(nodes_a) + length(nodes_b)) * 3;

        m_variables.resize(nb_dofs);

        auto variable_it = m_variables.begin();

        for (index i = 0; i < length(nodes_a); i++) {
            *(variable_it++) = nodes_a[i]->x();
            *(variable_it++) = nodes_a[i]->y();
            *(variable_it++) = nodes_a[i]->z();
        }
        
        for (index i = 0; i < length(nodes_b); i++) {
            *(variable_it++) = nodes_b[i]->x();
            *(variable_it++) = nodes_b[i]->y();
            *(variable_it++) = nodes_b[i]->z();
        }
    }

    template <int TOrder>
    double compute(Ref<Vector> g, Ref<Matrix> h) const
    {
        double f = 0;
        g.setZero();
        h.setZero();

        for (const auto& [shape_functions_a, shape_functions_b, weight] : m_data) {
            Vector3D point_a = ref_geometry(m_nodes_a, shape_functions_a, 0);
            Vector3D point_b = ref_geometry(m_nodes_b, shape_functions_b, 0);

            const Vector3D delta = point_b - point_a;

            f += delta.squaredNorm() * weight / 2;
            
            if constexpr(TOrder > 0) {
                index value_index = 0;

                for (index i = 0; i < length(m_nodes_a); i++) {
                    g(value_index++) += delta[0] * shape_functions_a(0, i) * weight;
                    g(value_index++) += delta[1] * shape_functions_a(0, i) * weight;
                    g(value_index++) += delta[2] * shape_functions_a(0, i) * weight;
                }

                for (index i = 0; i < length(m_nodes_b); i++) {
                    g(value_index++) -= delta[0] * shape_functions_b(0, i) * weight;
                    g(value_index++) -= delta[1] * shape_functions_b(0, i) * weight;
                    g(value_index++) -= delta[2] * shape_functions_b(0, i) * weight;
                }
            }

            if constexpr(TOrder > 1) {
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

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Holder = Pointer<Type>;
        using Base = Objective;

        py::class_<Type, Base, Holder>(m, "IgaDisplacementCoupling")
            .def(py::init<std::vector<Pointer<Node>>, std::vector<Pointer<Node>>, std::vector<Data>>())
        ;
    }
};

} // namespace eqlib