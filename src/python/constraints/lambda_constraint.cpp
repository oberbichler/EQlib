#include "../common.h"

#include <eqlib/variable.h>
#include <eqlib/constraints/lambda_constraint.h>

void bind_lambda_constraint(py::module &m, py::module &s)
{
    namespace py = pybind11;
    using namespace eqlib;
    using namespace eqlib::constraints;
    using namespace pybind11::literals;

    auto cls = bind_constraint<LambdaConstraint>(s, "LambdaConstraint");

    cls
        .def(py::init<Equations, Variables, LambdaConstraint::Fn>(), "equations"_a, "variables"_a, "fn"_a);

    m.def("constraint", [](const Equations& equations, const Variables& variables) {
        auto factory = [=](LambdaConstraint::Fn fn) {
            return LambdaConstraint(equations, variables, fn);
        };
        return py::cpp_function(factory, "fn"_a);
    }, "equations"_a, "variables"_a);
}