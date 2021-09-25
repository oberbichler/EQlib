#include "common.h"

#include <eqlib/problem.h>

void bind_request(py::module& m)
{
    using namespace eqlib;

    py::enum_<Request>(m, "Request", py::arithmetic())
        .value("F", Request::F)
        .value("DF", Request::Df)
        .value("G", Request::G)
        .value("DG", Request::Dg)
        .value("HF", Request::Hf)
        .value("HG", Request::Hg)
        .export_values();
}