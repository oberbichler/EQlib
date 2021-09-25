#include "common.h"

void bind_constraint(py::module_&);
void bind_problem(py::module_&);
void bind_equation(py::module_&);
void bind_node(py::module_&);
void bind_objective(py::module_&);
void bind_sparse_structure(py::module& m);
void bind_variable(py::module_&);

void bind_lambda_constraint(py::module_&, py::module_&);
void bind_lambda_objective(py::module_&, py::module_&);
void bind_spring(py::module_&);

PYBIND11_MODULE(eqlib_ext, m)
{
    m.doc() = "EQlib by Thomas Oberbichler";
    m.attr("__author__") = "Thomas Oberbichler";
    m.attr("__copyright__") = "Copyright (c) 2019-2021, Thomas Oberbichler";
    m.attr("__version__") = EQLIB_VERSION;
    m.attr("__email__") = "thomas.oberbichler@gmail.com";
    m.attr("__status__") = "Development";

    m.def("show_config", [] {
        MKLVersion Version;
 
        mkl_get_version(&Version);

        py::print("EQlib", EQLIB_VERSION, "with MKL", Version.MajorVersion, Version.MinorVersion, Version.UpdateVersion);
    });

    bind_constraint(m);
    bind_problem(m);
    bind_equation(m);
    bind_node(m);
    bind_objective(m);
    bind_sparse_structure(m);
    bind_variable(m);

    auto objectives_submodule = m.def_submodule("objectives");

    bind_lambda_constraint(m, objectives_submodule);
    bind_lambda_objective(m, objectives_submodule);
    bind_spring(objectives_submodule);
}