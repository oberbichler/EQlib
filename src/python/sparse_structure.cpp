#include "common.h"

#include <eqlib/sparse_structure.h>

void bind_sparse_structure(py::module_ &m)
{
    using namespace eqlib;

    auto csr_structure_cls = bind<CsrStructure>(m, "CsrStructure");

    csr_structure_cls
        .def(py::init<int, int, std::vector<int>, std::vector<int>>(), "rows"_a, "cols"_a, "ia"_a, "ja"_a)
        .def_property_readonly("rows", &CsrStructure::rows)
        .def_property_readonly("cols", &CsrStructure::cols)
        .def_property_readonly("nb_nonzeros", py::overload_cast<>(&CsrStructure::nb_nonzeros, py::const_))
        .def_property_readonly("ia", py::overload_cast<>(&CsrStructure::ia, py::const_))
        .def_property_readonly("ja", py::overload_cast<>(&CsrStructure::ja, py::const_))
        .def("get_index", &CsrStructure::get_index, "row"_a, "col"_a);
}