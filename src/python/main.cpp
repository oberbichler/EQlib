#include "common.h"

// eqlib::
void bind_constraint(py::module_&);
void bind_problem(py::module_&);
void bind_equation(py::module_&);
void bind_node(py::module_&);
void bind_objective(py::module_&);
void bind_parameter(py::module_&);
void bind_request(py::module_&);
void bind_sparse_structure(py::module_&);

// eqlib::objectives
void bind_lambda_objective(py::module_&, py::module_&);
void bind_spring(py::module_&);

// eqlib::constraints
void bind_lambda_constraint(py::module_&, py::module_&);

// eqlib::linear_solvers
void bind_linear_solver(py::module_&, py::module_&);
#ifdef EQLIB_USE_MKL
void bind_pardiso_ldlt(py::module_&, py::module_&);
#endif
void bind_simplicial_ldlt(py::module_&, py::module_&);

PYBIND11_MODULE(eqlib_ext, m)
{
    m.doc() = "EQlib by Thomas Oberbichler";
    m.attr("__author__") = "Thomas Oberbichler";
    m.attr("__copyright__") = "Copyright (c) 2019-2021, Thomas Oberbichler";
    m.attr("__version__") = EQLIB_VERSION;
    m.attr("__email__") = "thomas.oberbichler@gmail.com";
    m.attr("__status__") = "Development";

    // eqlib::

    m.def("show_config", [] {
        py::print("EQlib", EQLIB_VERSION);
    });

    m.attr("LO") = eqlib::LOWEST;
    m.attr("HI") = eqlib::HIGHEST;

    bind_request(m);
    bind_constraint(m);
    bind_problem(m);
    bind_equation(m);
    bind_node(m);
    bind_objective(m);
    bind_parameter(m);
    bind_sparse_structure(m);

    // eqlib::objectives
    auto objectives_submodule = m.def_submodule("objectives");

    bind_lambda_objective(m, objectives_submodule);
    bind_spring(objectives_submodule);

    // eqlib::constraints
    auto constraints_submodule = m.def_submodule("constraints");
    bind_lambda_constraint(m, constraints_submodule);

    // eqlib::linear_solvers
    auto linear_solvers_submodule = m.def_submodule("linear_solvers");
    bind_linear_solver(m, linear_solvers_submodule);
#ifdef EQLIB_USE_MKL
    bind_pardiso_ldlt(m, linear_solvers_submodule);
#endif
    bind_simplicial_ldlt(m, linear_solvers_submodule);
}