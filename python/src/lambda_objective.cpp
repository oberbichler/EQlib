#include "common.h"

#include <eqlib/variable.h>
#include <eqlib/lambda_objective.h>

void bind_lambda_objective(py::module &m, py::module &s)
{
    namespace py = pybind11;
    using namespace eqlib;
    using namespace pybind11::literals;

    auto cls = bind_objective<LambdaObjective>(s, "LambdaObjective");

    cls
        .def(py::init<Variables, LambdaObjective::Fn>(), "variables"_a, "fn"_a);

    m.def("objective", [](const Variables& variables) {
        auto factory = [=](LambdaObjective::Fn fn) {
            return LambdaObjective(variables, fn);
        };
        return py::cpp_function(factory, "fn"_a);
    }, "variables"_a);
}