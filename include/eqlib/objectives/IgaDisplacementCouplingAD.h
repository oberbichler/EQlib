#pragma once

#include <Eigen/Geometry>

#include <eqlib/Objective.h>
#include <eqlib/Node.h>
#include <eqlib/Variable.h>

#include <hyperjet/hyperjet.h>

#include <vector>

namespace eqlib {

class IgaDisplacementCouplingAD : public Objective
{
private:    // types
    using Type = IgaDisplacementCouplingAD;
    using Jet = hyperjet::Jet<double>;
    using HyperJet = hyperjet::HyperJet<double>;
    using Jet3D = Eigen::Matrix<Jet, 3, 1>;
    using HyperJet3D = Eigen::Matrix<HyperJet, 3, 1>;

    using Data = std::tuple<Matrix, Matrix, double>;

private:    // variables
    std::vector<Pointer<Node>> m_nodes_a;
    std::vector<Pointer<Node>> m_nodes_b;
    std::vector<Data> m_data;
    double m_weight;

    template <int TOrder>
    auto ref_geometry(const std::vector<Pointer<Node>>& nodes, Ref<const Matrix> shape_functions, const index i) const
    {
        Vector3D value = Vector3D::Zero();

        for (index j = 0; j < length(nodes); j++) {
            value += nodes[j]->ref_location() * shape_functions(i, j);
        }

        if constexpr(TOrder == 0) {
            return value;
        }

        if constexpr(TOrder == 1) {
            Jet3D jet;

            for (index k = 0; k < 3; k++) {
                jet[k] = Jet(value(k), length(nodes) * 3);
            }

            for (index j = 0; j < length(nodes); j++) {
                jet(0).g(j * 3 + 0) = shape_functions(i, j);
                jet(1).g(j * 3 + 1) = shape_functions(i, j);
                jet(2).g(j * 3 + 2) = shape_functions(i, j);
            }

            return jet;
        }

        if constexpr(TOrder == 2) {
            HyperJet3D hyper_jet;

            for (index k = 0; k < 3; k++) {
                hyper_jet[k] = HyperJet(value(k), length(nodes) * 3);
            }

            for (index j = 0; j < length(nodes); j++) {
                hyper_jet(0).g(j * 3 + 0) = shape_functions(i, j);
                hyper_jet(1).g(j * 3 + 1) = shape_functions(i, j);
                hyper_jet(2).g(j * 3 + 2) = shape_functions(i, j);
            }

            return hyper_jet;
        }
    }

    template <int TOrder>
    auto act_geometry(const std::vector<Pointer<Node>>& nodes, Ref<const Matrix> shape_functions, const index i) const
    {
        Vector3D value = Vector3D::Zero();

        for (index j = 0; j < length(nodes); j++) {
            value += nodes[j]->act_location() * shape_functions(i, j);
        }

        if constexpr(TOrder == 0) {
            return value;
        }

        if constexpr(TOrder == 1) {
            Jet3D jet;

            for (index k = 0; k < 3; k++) {
                jet[k] = Jet(value(k), length(nodes) * 3);
            }

            for (index j = 0; j < length(nodes); j++) {
                jet(0).g(j * 3 + 0) = shape_functions(i, j);
                jet(1).g(j * 3 + 1) = shape_functions(i, j);
                jet(2).g(j * 3 + 2) = shape_functions(i, j);
            }

            return jet;
        }

        if constexpr(TOrder == 2) {
            HyperJet3D hyper_jet;

            for (index k = 0; k < 3; k++) {
                hyper_jet[k] = HyperJet(value(k), length(nodes) * 3);
            }

            for (index j = 0; j < length(nodes); j++) {
                hyper_jet(0).g(j * 3 + 0) = shape_functions(i, j);
                hyper_jet(1).g(j * 3 + 1) = shape_functions(i, j);
                hyper_jet(2).g(j * 3 + 2) = shape_functions(i, j);
            }

            return hyper_jet;
        }
    }

public:     // constructor
    IgaDisplacementCouplingAD(
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
            auto point_a = ref_geometry<TOrder>(m_nodes_a, shape_functions_a, 0);
            auto point_b = ref_geometry<TOrder>(m_nodes_b, shape_functions_b, 0);

            if constexpr(TOrder > 0) {
                const index nb_dofs_a = length(m_nodes_a) * 3;
                const index nb_dofs_b = length(m_nodes_b) * 3;

                for (index i = 0; i < 3; i++) {
                    point_a(i) = point_a(i).enlarge(0, nb_dofs_b);
                    point_b(i) = point_b(i).enlarge(nb_dofs_a, 0);
                }
            }

            const auto ab = point_b - point_a;

            const auto p = ab.dot(ab) * m_weight / 2;

            f += hyperjet::explode_add(p, g, h);
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

        py::class_<Type, Base, Holder>(m, "IgaDisplacementCouplingAD")
            .def(py::init<std::vector<Pointer<Node>>, std::vector<Pointer<Node>>, std::vector<Data>>())
        ;
    }
};

} // namespace eqlib