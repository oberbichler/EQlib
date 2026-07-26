#include "bindings/bindings.h"

#include <eqlib/Info.h>

PYBIND11_MODULE(eqlib, m)
{
    m.doc() = "eqlib by Thomas Oberbichler";
    m.attr("__author__") = "Thomas Oberbichler";
    m.attr("__copyright__") = "Copyright (c) 2018-2026, Thomas Oberbichler";
    m.attr("__version__") = eqlib::Info::version();
    m.attr("__email__") = "thomas.oberbichler@gmail.com";
    m.attr("__status__") = "Development";

    m.attr("_GIT_COMMIT_HASH") = eqlib::Info::git_commit_hash();
    m.attr("_USE_BLAS") = eqlib::Info::use_blas();
    m.attr("_USE_MKL") = eqlib::Info::use_mkl();

    // --- core

    eqlib::register_armijo(m);
    eqlib::register_equation(m);
    eqlib::register_variable(m);
    eqlib::register_parameter(m);
    eqlib::register_constraint(m);
    eqlib::register_objective(m);
    eqlib::register_lambda_constraint(m);
    eqlib::register_lambda_objective(m);
    eqlib::register_problem(m);
    eqlib::register_log(m);
    eqlib::register_node(m);
    eqlib::register_timer(m);
    eqlib::register_newton_raphson(m);
    eqlib::register_steepest_descent(m);
    eqlib::register_linear_solver(m);
    eqlib::register_simplicial_ldlt(m);
    eqlib::register_sparse_lu(m);

#ifdef EQLIB_USE_MKL
    eqlib::register_pardiso_ldlt(m);
#endif

    eqlib::register_sparse_structures(m);

    // --- objectives

    eqlib::register_iga_normal_distance_ad(m);
    eqlib::register_iga_membrane_3p_ad(m);
    eqlib::register_iga_point_distance(m);
    eqlib::register_iga_point_distance_ad(m);
    eqlib::register_iga_point_location(m);
    eqlib::register_iga_rotation_coupling_ad(m);
    eqlib::register_iga_shell_3p_ad(m);
}
