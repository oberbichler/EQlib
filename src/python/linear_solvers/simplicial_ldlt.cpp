#include "../common.h"

#include <eqlib/linear_solvers/simplicial_ldlt.h>

void bind_simplicial_ldlt(py::module_ &m, py::module_ &s)
{
    namespace py = pybind11;
    using namespace eqlib;
    using namespace eqlib::linear_solvers;
    using namespace pybind11::literals;

    py::class_<SimplicialLDLT, LinearSolver, Pointer<SimplicialLDLT>>(s, "SimplicialLDLT")
        .def(py::init<>());
}