#pragma once

#include "../Node.h"
#include "../Objective.h"

#include <hyperjet/hyperjet.h>

namespace eqlib {

namespace iga_utilities {

Vector3D evaluate_ref_geometry(const std::vector<Pointer<Node>>& nodes, Ref<const Vector> shape_functions)
{
    Vector3D value = Vector3D::Zero();

    for (index i = 0; i < length(nodes); i++) {
        value += nodes[i]->ref_location() * shape_functions(i);
    }

    return value;
}

Vector3D evaluate_act_geometry(const std::vector<Pointer<Node>>& nodes, Ref<const Vector> shape_functions)
{
    Vector3D value = Vector3D::Zero();

    for (index i = 0; i < length(nodes); i++) {
        value += nodes[i]->act_location() * shape_functions(i);
    }

    return value;
}

auto evaluate_act_geometry_hj(const std::vector<Pointer<Node>>& nodes, Ref<const Vector> shape_functions)
{
    using Space = hyperjet::Space<2, double, -1>;

    Vector3D xyz = evaluate_act_geometry(nodes, shape_functions);

    Space::Vector<3> value;

    const index nb_variables = length(nodes) * 3;

    value[0] = Space::Scalar::constant(nb_variables, xyz[0]);
    value[1] = Space::Scalar::constant(nb_variables, xyz[1]);
    value[2] = Space::Scalar::constant(nb_variables, xyz[2]);

    for (index i = 0; i < length(shape_functions); i++) {
        Space::set_g(value[0], i * 3 + 0, shape_functions[i]);
        Space::set_g(value[1], i * 3 + 1, shape_functions[i]);
        Space::set_g(value[2], i * 3 + 2, shape_functions[i]);
    }

    return value;
}

auto evaluate_act_geometry_hj_a(const std::vector<Pointer<Node>>& nodes, Ref<const Vector> shape_functions, index nb_variables)
{
    using Space = hyperjet::Space<2, double, -1>;

    Vector3D xyz = evaluate_act_geometry(nodes, shape_functions);

    Space::Vector<3> value;

    value[0] = Space::Scalar::constant(nb_variables, xyz[0]);
    value[1] = Space::Scalar::constant(nb_variables, xyz[1]);
    value[2] = Space::Scalar::constant(nb_variables, xyz[2]);

    for (index i = 0; i < length(shape_functions); i++) {
        Space::set_g(value[0], i * 3 + 0, shape_functions[i]);
        Space::set_g(value[1], i * 3 + 1, shape_functions[i]);
        Space::set_g(value[2], i * 3 + 2, shape_functions[i]);
    }

    return value;
}

auto evaluate_act_geometry_hj_b(const std::vector<Pointer<Node>>& nodes, Ref<const Vector> shape_functions, index nb_variables)
{
    using Space = hyperjet::Space<2, double, -1>;

    Vector3D xyz = evaluate_act_geometry(nodes, shape_functions);

    Space::Vector<3> value;

    value[0] = Space::Scalar::constant(nb_variables, xyz[0]);
    value[1] = Space::Scalar::constant(nb_variables, xyz[1]);
    value[2] = Space::Scalar::constant(nb_variables, xyz[2]);

    const index offset = nb_variables - length(nodes) * 3;

    for (index i = 0; i < length(shape_functions); i++) {
        Space::set_g(value[0], offset + i * 3 + 0, shape_functions[i]);
        Space::set_g(value[1], offset + i * 3 + 1, shape_functions[i]);
        Space::set_g(value[2], offset + i * 3 + 2, shape_functions[i]);
    }

    return value;
}

} // namespace iga_utilities

} // namespace eqlib