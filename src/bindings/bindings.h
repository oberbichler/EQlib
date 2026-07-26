#pragma once

#include <pybind11/pybind11.h>

namespace eqlib {

void register_armijo(pybind11::module_& m);
void register_equation(pybind11::module_& m);
void register_variable(pybind11::module_& m);
void register_parameter(pybind11::module_& m);
void register_constraint(pybind11::module_& m);
void register_objective(pybind11::module_& m);
void register_lambda_constraint(pybind11::module_& m);
void register_lambda_objective(pybind11::module_& m);
void register_problem(pybind11::module_& m);
void register_log(pybind11::module_& m);
void register_node(pybind11::module_& m);
void register_timer(pybind11::module_& m);
void register_newton_raphson(pybind11::module_& m);
void register_steepest_descent(pybind11::module_& m);
void register_linear_solver(pybind11::module_& m);
void register_simplicial_ldlt(pybind11::module_& m);
void register_sparse_lu(pybind11::module_& m);
#ifdef EQLIB_USE_MKL
void register_pardiso_ldlt(pybind11::module_& m);
#endif
void register_sparse_structures(pybind11::module_& m);
void register_iga_normal_distance_ad(pybind11::module_& m);
void register_iga_membrane_3p_ad(pybind11::module_& m);
void register_iga_point_distance(pybind11::module_& m);
void register_iga_point_distance_ad(pybind11::module_& m);
void register_iga_point_location(pybind11::module_& m);
void register_iga_rotation_coupling_ad(pybind11::module_& m);
void register_iga_shell_3p_ad(pybind11::module_& m);

} // namespace eqlib
