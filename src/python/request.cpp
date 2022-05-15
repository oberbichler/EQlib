#include "common.h"

#include <eqlib/problem.h>

void bind_request(py::module_& m)
{
    using namespace eqlib;

    py::enum_<Request>(m, "Request", py::arithmetic())
        .value("F", Request::F)
        .value("DF", Request::DF)
        .value("HF", Request::HF)
        .value("G", Request::G)
        .value("DG", Request::DG)
        .value("HG", Request::HG)
        .value("HM", Request::HM)
        .value("All", Request::All)
        .export_values();
}