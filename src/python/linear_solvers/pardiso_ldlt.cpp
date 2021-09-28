#ifdef EQLIB_USE_MKL

#include "../common.h"

#include <eqlib/linear_solvers/pardiso_ldlt.h>

void bind_pardiso_ldlt(py::module &m, py::module &s)
{
    namespace py = pybind11;
    using namespace eqlib;
    using namespace eqlib::linear_solvers;
    using namespace pybind11::literals;

    py::class_<PardisoLDLT, LinearSolver, Pointer<PardisoLDLT>>(s, "PardisoLDLT")
        .def(py::init<>());
}

#endif