#include "../common.h"

#include <eqlib/parameter.h>
#include <eqlib/constraints/lambda_constraint.h>

void bind_lambda_constraint(py::module_ &m, py::module_ &s)
{
    namespace py = pybind11;
    using namespace eqlib;
    using namespace eqlib::constraints;
    using namespace pybind11::literals;

    auto cls = bind_constraint<LambdaConstraint>(s, "LambdaConstraint");

    cls
        .def(py::init<Equations, Parameters, LambdaConstraint::Fn>(), "equations"_a, "parameters"_a, "fn"_a);

    m.def("constraint", [](const Equations& equations, const Parameters& parameters) {
        auto factory = [=](LambdaConstraint::Fn fn) {
            return LambdaConstraint(equations, parameters, fn);
        };
        return py::cpp_function(factory, "fn"_a);
    }, "equations"_a, "parameters"_a);
}