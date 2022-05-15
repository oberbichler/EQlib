#include "../common.h"

#include <eqlib/parameter.h>
#include <eqlib/objectives/lambda_objective.h>

void bind_lambda_objective(py::module_ &m, py::module_ &s)
{
    namespace py = pybind11;
    using namespace eqlib;
    using namespace eqlib::objectives;
    using namespace pybind11::literals;

    auto cls = bind_objective<LambdaObjective>(s, "LambdaObjective");

    cls
        .def(py::init<Parameters, LambdaObjective::Fn>(), "parameters"_a, "fn"_a);

    m.def("objective", [](const Parameters& parameters) {
        auto factory = [=](LambdaObjective::Fn fn) {
            return LambdaObjective(parameters, fn);
        };
        return py::cpp_function(factory, "fn"_a);
    }, "parameters"_a);
}