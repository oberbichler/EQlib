#include "common.h"

#include <eqlib/objective.h>

void bind_objective(py::module_ &m)
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

    bind<Objective, PyObjective>(m, "Objective")
        .def(py::init<>())
        .def_property("is_active", &Objective::is_active, &Objective::set_is_active)
        .def_property("name", &Objective::name, &Objective::set_name)
        .def_property_readonly("nb_parameters", &Objective::nb_parameters)
        .def_property("parameters", &Objective::parameters, &Objective::set_parameters)
        .def_property("parameter_values", &Objective::parameter_values, &Objective::set_parameter_values);
}