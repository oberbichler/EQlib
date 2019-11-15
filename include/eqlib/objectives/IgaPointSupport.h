#pragma once

#include <eqlib/Node.h>
#include <eqlib/Objective.h>
#include <eqlib/Variable.h>

#include <vector>

namespace eqlib {

class IgaPointSupport : public Objective
{
private:    // types

private:    // variables
    std::vector<Pointer<Node>> m_nodes;
    std::vector<std::tuple<Matrix, Vector3D>> m_data;
    double m_weight;

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
        std::vector<std::tuple<Matrix, Vector3D>> shape_functions,
        double weight)
    : m_data(shape_functions)
    , m_nodes(nodes)
    , m_weight(weight)
    {
        m_variables.resize(length(nodes) * 3);
        for (index i = 0; i < length(nodes); i++) {
            m_variables[i * 3 + 0] = nodes[i]->x();
            m_variables[i * 3 + 1] = nodes[i]->y();
            m_variables[i * 3 + 2] = nodes[i]->z();
        }
    }

    double compute(Ref<Vector> g, Ref<Matrix> h) const override
    {
        double f = 0;
            g.setZero();
            h.setZero();

        for (const auto& [shape_functions, target] : m_data) {
            const Vector3D act_x = act_geometry(0, shape_functions);

            const Vector3D delta = target - act_x;

            f += delta.dot(delta) * m_weight / 2;

            for (index i = 0; i < length(m_nodes); i++) {
                g(i * 3 + 0) += delta[0] * shape_functions(0, i) * m_weight;
                g(i * 3 + 1) += delta[1] * shape_functions(0, i) * m_weight;
                g(i * 3 + 2) += delta[2] * shape_functions(0, i) * m_weight;
            }

            for (index i = 0; i < length(m_nodes); i++) {
                for (index j = 0; j < length(m_nodes); j++) {
                    h(i * 3 + 0, j * 3 + 0) += shape_functions(0, i) * shape_functions(0, j) * m_weight;
                    h(i * 3 + 1, j * 3 + 1) += shape_functions(0, i) * shape_functions(0, j) * m_weight;
                    h(i * 3 + 2, j * 3 + 2) += shape_functions(0, i) * shape_functions(0, j) * m_weight;
                }
            }
        }
        
        return f;
    }

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Type = IgaPointSupport;
        using Holder = Pointer<Type>;
        using Base = Objective;

        py::class_<Type, Base, Holder>(m, "IgaPointSupport")
            .def(py::init<std::vector<Pointer<Node>>, std::vector<std::tuple<Matrix, Vector3D>>, double>())
        ;
    }
};

} // namespace eqlib