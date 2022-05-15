#include "common.h"

#include <eqlib/problem.h>

void bind_problem(py::module_& m)
{
    using namespace eqlib;

    py::object scipy_sparse = py::module_::import("scipy.sparse");
    py::object csr_matrix = scipy_sparse.attr("csr_matrix");

    bind<Problem>(m, "Problem")
        .def(py::init<>(&Problem::create), "objectives"_a = py::list(), "constraints"_a = py::list(), py::keep_alive<1, 2>(), py::keep_alive<1, 3>())
        .def_static("like", &Problem::like, "prototype"_a, "objectives"_a = py::list(), "constraints"_a = py::list(), py::keep_alive<1, 2>(), py::keep_alive<1, 3>())
        .def("compute_to", py::overload_cast<Ref<Vector>, Ref<Vector>, Ref<Vector>, Ref<Vector>>(&Problem::compute), "g"_a = py::array(), "df"_a = py::array(), "dg"_a = py::array(), "hm"_a = py::array(), py::call_guard<py::gil_scoped_release>())
        .def("compute", py::overload_cast<Request>(&Problem::compute), "request"_a = Request::All, py::call_guard<py::gil_scoped_release>())
        .def_property_readonly("structure_dg", &Problem::structure_dg)
        .def_property_readonly("structure_hm", &Problem::structure_hm)
        .def_property_readonly("nb_objectives", &Problem::nb_objectives)
        .def_property_readonly("nb_constraints", &Problem::nb_constraints)
        .def_property_readonly("nb_equations", &Problem::nb_equations)
        .def_property_readonly("equations", &Problem::equations)
        .def_property_readonly("nb_variables", &Problem::nb_variables)
        .def_property_readonly("variables", &Problem::variables)
        .def_property_readonly("nb_constants", &Problem::nb_constants)
        .def_property_readonly("constants", &Problem::constants)
        .def_property_readonly("buffer_size", &Problem::buffer_size)
        .def_property_readonly("nb_nonzeros_dg", &Problem::nb_nonzeros_dg)
        .def_property_readonly("nb_nonzeros_hm", &Problem::nb_nonzeros_hm)
        .def("equation_index", &Problem::equation_index, "equation"_a)
        .def("variable_index", &Problem::variable_index, "variable"_a)
        .def("constant_index", &Problem::constant_index, "constant"_a)
        .def_property("x", &Problem::x, &Problem::set_x)
        .def_property("x_lower_bounds", &Problem::x_lower_bounds, &Problem::set_x_lower_bounds)
        .def_property("x_upper_bounds", &Problem::x_upper_bounds, &Problem::set_x_upper_bounds)
        .def("add_x", &Problem::add_x, "value"_a)
        .def("sub_x", &Problem::sub_x, "value"_a)
        .def("eval", &Problem::eval, "x"_a)
        .def_property_readonly("linear_solver", &Problem::linear_solver)
        .def_property_readonly("f", &Problem::f)
        .def_property_readonly("g", &Problem::g)
        .def_property_readonly("df", &Problem::df)
        .def_property_readonly("dg", [=](Problem& self) {
            auto data = std::make_tuple(self.dg_values(), self.structure_dg().ja(), self.structure_dg().ia());
            auto shape = std::make_pair(self.nb_equations(), self.nb_variables());
            return csr_matrix(data, shape).release();
        })
        .def_property_readonly("hm", [=](Problem& self) {
            auto data = std::make_tuple(self.hm_values(), self.structure_hm().ja(), self.structure_hm().ia());
            auto shape = std::make_pair(self.nb_variables(), self.nb_variables());
            return csr_matrix(data, shape).release();
        })
        .def_property("parallel_tasks", &Problem::parallel_tasks, &Problem::set_parallel_tasks)
        .def_property_readonly("dg_values", &Problem::dg_values)
        .def_property_readonly("hm_values", &Problem::hm_values)
        .def_property("hm_diagonal", &Problem::hm_diagonal, &Problem::set_hm_diagonal)
        .def("hm_add_diagonal", &Problem::hm_add_diagonal)
        .def("newton_step", &Problem::newton_step)
        .def("hm_inv_v", &Problem::hm_inv_v, "v"_a)
        .def("hm_v", &Problem::hm_v, "v"_a)
        .def("scale_hm", &Problem::scale_hm, "factor"_a)
        .def_property_readonly("hm_norm_inf", &Problem::hm_norm_inf);
}