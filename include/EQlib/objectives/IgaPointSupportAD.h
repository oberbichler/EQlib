#pragma once

#include <Eigen/Geometry>

#include <EQlib/Objective.h>
#include <EQlib/Node.h>
#include <EQlib/Variable.h>

#include <hyperjet/hyperjet.h>

#include <vector>

namespace EQlib {

class IgaPointSupportAD : public Objective
{
private:    // types
    using Jet = hyperjet::Jet<double>;
    using HyperJet = hyperjet::HyperJet<double>;

    using Jet3D = Eigen::Matrix<Jet, 3, 1>;
    using HyperJet3D = Eigen::Matrix<HyperJet, 3, 1>;

private:    // variables
    std::vector<Pointer<Node>> m_nodes;
    Matrix m_shape_functions;
    double m_weight;
    Vector3D m_target;

    template <int TOrder>
    auto ref_geometry(index i) const
    {
        Vector3D value = Vector3D::Zero();

        for (index j = 0; j < length(m_nodes); j++) {
            value += m_nodes[j]->ref_location() * m_shape_functions(i, j);
        }

        if constexpr(TOrder == 0) {
            return value;
        }

        if constexpr(TOrder == 1) {
            Jet3D jet;

            for (index k = 0; k < 3; k++) {
                jet[k] = Jet(value(k), length(m_nodes) * 3);
            }

            for (index j = 0; j < length(m_nodes); j++) {
                jet(0).g(j * 3 + 0) = m_shape_functions(i, j);
                jet(1).g(j * 3 + 1) = m_shape_functions(i, j);
                jet(2).g(j * 3 + 2) = m_shape_functions(i, j);
            }

            return jet;
        }

        if constexpr(TOrder == 2) {
            HyperJet3D hyper_jet;

            for (index k = 0; k < 3; k++) {
                hyper_jet[k] = HyperJet(value(k), length(m_nodes) * 3);
            }

            for (index j = 0; j < length(m_nodes); j++) {
                hyper_jet(0).g(j * 3 + 0) = m_shape_functions(i, j);
                hyper_jet(1).g(j * 3 + 1) = m_shape_functions(i, j);
                hyper_jet(2).g(j * 3 + 2) = m_shape_functions(i, j);
            }

            return hyper_jet;
        }
    }

    template <int TOrder>
    auto act_geometry(index i) const
    {
        Vector3D value = Vector3D::Zero();

        for (index j = 0; j < length(m_nodes); j++) {
            value += m_nodes[j]->act_location() * m_shape_functions(i, j);
        }

        if constexpr(TOrder == 0) {
            return value;
        }

        if constexpr(TOrder == 1) {
            Jet3D jet;

            for (index k = 0; k < 3; k++) {
                jet[k] = Jet(value(k), length(m_nodes) * 3);
            }

            for (index j = 0; j < length(m_nodes); j++) {
                jet(0).g(j * 3 + 0) = m_shape_functions(i, j);
                jet(1).g(j * 3 + 1) = m_shape_functions(i, j);
                jet(2).g(j * 3 + 2) = m_shape_functions(i, j);
            }

            return jet;
        }

        if constexpr(TOrder == 2) {
            HyperJet3D hyper_jet;

            for (index k = 0; k < 3; k++) {
                hyper_jet[k] = HyperJet(value(k), length(m_nodes) * 3);
            }

            for (index j = 0; j < length(m_nodes); j++) {
                hyper_jet(0).g(j * 3 + 0) = m_shape_functions(i, j);
                hyper_jet(1).g(j * 3 + 1) = m_shape_functions(i, j);
                hyper_jet(2).g(j * 3 + 2) = m_shape_functions(i, j);
            }

            return hyper_jet;
        }
    }

public:     // constructor
    IgaPointSupportAD(
        std::vector<Pointer<Node>> nodes,
        Matrix shape_functions,
        Vector3D target,
        double weight)
    : m_shape_functions(shape_functions)
    , m_nodes(nodes)
    , m_target(target)
    , m_weight(weight)
    {
        m_variables.resize(length(nodes) * 3);
        for (index i = 0; i < length(nodes); i++) {
            m_variables[i * 3 + 0] = nodes[i]->x();
            m_variables[i * 3 + 1] = nodes[i]->y();
            m_variables[i * 3 + 2] = nodes[i]->z();
        }
    }

    template <int TOrder>
    double compute_(Ref<Vector> g, Ref<Matrix> h) const
    {
        // reference configuration

        // const auto ref_x = ref_geometry<0>(0);

        // actual configuration

        const auto act_x = act_geometry<TOrder>(0);

        // functional

        if constexpr(TOrder == 0) {
            const Vector3D u = m_target - act_x;

            const double p = u.dot(u) * 0.5 * m_weight;

            return hyperjet::explode(p, g, h);
        }
        if constexpr(TOrder == 1) {
            const Jet3D u = m_target - act_x;

            const Jet p = u.dot(u) * 0.5 * m_weight;

            return hyperjet::explode(p, g, h);
        }
        if constexpr(TOrder == 2) {
            const HyperJet3D u = m_target - act_x;

            const HyperJet p = u.dot(u) * 0.5 * m_weight;

            return hyperjet::explode(p, g, h);
        }
    }

    double compute(Ref<Vector> g, Ref<Matrix> h) const override
    {
        if (g.size() == 0) {
            return compute_<0>(g, h);
        } else if (h.size() == 0) {
            return compute_<1>(g, h);
        } else {
            return compute_<2>(g, h);
        }
    }

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Type = IgaPointSupportAD;
        using Holder = Pointer<Type>;
        using Base = Objective;

        py::class_<Type, Base, Holder>(m, "IgaPointSupportAD")
            .def(py::init<std::vector<Pointer<Node>>, Matrix, Vector3D, double>())
        ;
    }
};

} // namespace EQlib