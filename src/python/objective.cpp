#include "common.h"

#include <eqlib/objective.h>

void bind_objective(py::module &m)
{
    namespace py = pybind11;
    using namespace eqlib;
    using namespace pybind11::literals;

    struct PyObjective : public Objective
    {
        using Objective::Objective;

        double compute(Ref<Vector> df, Ref<Matrix> hm, const Request request) override
        {
            pybind11::gil_scoped_acquire acquire;
            PYBIND11_OVERLOAD_PURE(double, Objective, compute, df, hm, request);
        }
    };

    auto cls = bind<Objective, PyObjective>(m, "Objective");

    cls
        .def(py::init<>())
        .def_property("is_active", &Objective::is_active, &Objective::set_is_active)
        .def_property("name", &Objective::name, &Objective::set_name)
        .def_property_readonly("nb_variables", &Objective::nb_variables)
        .def_property("variables", &Objective::variables, &Objective::set_variables)
        .def_property("variable_values", &Objective::variable_values, &Objective::set_variable_values);
}