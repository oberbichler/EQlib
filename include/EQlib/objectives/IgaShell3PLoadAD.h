#pragma once

#include <Eigen/Geometry>

#include <EQlib/Objective.h>
#include <EQlib/Point.h>
#include <EQlib/Variable.h>

#include <hyperjet/hyperjet.h>

#include <vector>

namespace EQlib {

class IgaShell3PLoadAD : public Objective
{
private:    // types
    using Jet = hyperjet::Jet<double>;
    using HyperJet = hyperjet::HyperJet<double>;

    using Jet3D = Eigen::Matrix<Jet, 3, 1>;
    using HyperJet3D = Eigen::Matrix<HyperJet, 3, 1>;

private:    // variables
    std::vector<Pointer<Point>> m_nodes;
    Matrix m_shape_functions;
    Vector3D m_load;
    Vector3D m_x_ref;
    double m_weight;

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
    IgaShell3PLoadAD(
        std::vector<Pointer<Point>> nodes,
        Matrix shape_functions,
        Vector3D load,
        double weight)
    : m_shape_functions(shape_functions)
    , m_nodes(nodes)
    , m_load(load)
    , m_weight(weight)
    {
        m_variables.resize(length(nodes) * 3);
        for (index i = 0; i < length(nodes); i++) {
            m_variables[i * 3 + 0] = nodes[i]->x();
            m_variables[i * 3 + 1] = nodes[i]->y();
            m_variables[i * 3 + 2] = nodes[i]->z();
        }

        // reference configuration

        m_x_ref = act_geometry<0>(0);
    }

    double compute(Ref<Vector> g, Ref<Matrix> h) const override
    {
        // reference configuration

        const Vector3D ref_x = ref_geometry<0>(0);
        const Vector3D ref_a1 = ref_geometry<0>(1);
        const Vector3D ref_a2 = ref_geometry<0>(2);
        const double ref_da = ref_a1.cross(ref_a2).norm();

        // actual configuration

        const HyperJet3D act_x = act_geometry<2>(0);
        const HyperJet3D act_a1 = act_geometry<2>(1);
        const HyperJet3D act_a2 = act_geometry<2>(2);
        const HyperJet act_da = act_a1.cross(act_a2).norm();

        // functional

        const HyperJet3D u = act_x - ref_x;

        const HyperJet sum = u(0) * m_load(0) + u(1) * m_load(1) + u(2) * m_load(2);

        const HyperJet p = -sum * m_weight * ref_da;

        return hyperjet::explode(p, g, h);
    }

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Type = IgaShell3PLoadAD;
        using Holder = Pointer<Type>;
        using Base = Objective;

        py::class_<Type, Base, Holder>(m, "IgaShell3PLoadAD")
            .def(py::init<std::vector<Pointer<Point>>, Matrix, Vector3D, double>())
        ;
    }
};

} // namespace EQlib