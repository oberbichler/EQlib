#include "common.h"

#include <eqlib/problem.h>

void bind_problem(py::module& m)
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

    auto problem_cls = bind<Problem>(m, "Problem");

    py::object scipy_sparse = py::module::import("scipy.sparse");
    py::object csr_matrix = scipy_sparse.attr("csr_matrix");

    problem_cls
        .def(py::init<>(&Problem::create), "objectives"_a = py::list(), "constraints"_a = py::list(), py::keep_alive<1, 2>(), py::keep_alive<1, 3>())
        .def_static("like", &Problem::like, "prototype"_a, "objectives"_a = py::list(), "constraints"_a = py::list(), py::keep_alive<1, 2>(), py::keep_alive<1, 3>())
        .def("compute", py::overload_cast<Ref<Vector>, Ref<Vector>, Ref<Vector>, Ref<Vector>>(&Problem::compute), "g"_a = py::array(), "df"_a = py::array(), "dg"_a = py::array(), "hm"_a = py::array(), py::call_guard<py::gil_scoped_release>())
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
        .def_property_readonly("f", &Problem::f)
        .def_property_readonly("g", &Problem::g)
        .def_property_readonly("df", &Problem::df)
        .def_property_readonly("dg", [=](Problem& self) {
            auto data = std::make_tuple(self.dg(), self.structure_dg().ja(), self.structure_dg().ia());
            auto shape = std::make_pair(self.nb_equations(), self.nb_variables());
            return csr_matrix(data, shape).release();
        })
        .def_property_readonly("hm", [=](Problem& self) {
            auto data = std::make_tuple(self.hm(), self.structure_hm().ja(), self.structure_hm().ia());
            auto shape = std::make_pair(self.nb_variables(), self.nb_variables());
            return csr_matrix(data, shape).release();
        })
        .def_property("parallel_level", &Problem::parallel_level, &Problem::set_parallel_level);
}