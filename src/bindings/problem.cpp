#include "bindings.h"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <eqlib/Problem.h>

namespace eqlib {

void register_problem(pybind11::module_& m)
{
    using Type = Problem;
    namespace py = pybind11;
    using namespace pybind11::literals;

    using Holder = Pointer<Type>;

    const std::string name = "Problem";

    py::object scipy_sparse = py::module::import("scipy.sparse");
    py::object csr_matrix = scipy_sparse.attr("csr_matrix");

    py::class_<Type, Holder>(m, name.c_str())
        // constructors
        .def(py::init<Type::ElementsF, Type::ElementsG, int, int>(), "objective"_a = py::list(), "constraints"_a = py::list(),
            "nb_threads"_a = 1, "grainsize"_a = 100, py::keep_alive<1, 2>(), py::keep_alive<1, 3>())
        // read-only properties
        .def_property_readonly("is_constrained", &Type::is_constrained)
        .def_property_readonly("equations", &Type::equations)
        .def_property_readonly("variables", &Type::variables)
        .def_property_readonly("g", py::overload_cast<>(&Type::g))
        .def_property_readonly("df", py::overload_cast<>(&Type::df))
        .def_property_readonly("dg", [=](Type& self) {
            return csr_matrix(
                std::make_tuple(self.dg_values(), self.dg_indices(), self.dg_indptr()),
                std::make_pair(self.nb_equations(), self.nb_variables()))
                .release();
        })
        .def_property_readonly("structure_dg", &Type::structure_dg)
        .def_property_readonly("dg_values", py::overload_cast<>(&Type::dg_values))
        .def_property_readonly("dg_indptr", &Type::dg_indptr)
        .def_property_readonly("dg_indices", &Type::dg_indices)
        .def_property_readonly("hm", [=](Type& self) {
            return csr_matrix(
                std::make_tuple(self.hm_values(), self.hm_indices(), self.hm_indptr()),
                std::make_pair(self.nb_variables(), self.nb_variables()))
                .release();
        })
        .def_property_readonly("general_hm", [=](Type& self) {
            const auto [structure, values] = self.structure_hm().to_general(self.hm_values());
            return csr_matrix(
                std::make_tuple(values, structure.ja(), structure.ia()),
                std::make_pair(self.nb_variables(), self.nb_variables()),
                "copy"_a=true)
                .release();
        })
        .def_property_readonly("structure_hm", &Type::structure_hm)
        .def_property_readonly("hm_values", py::overload_cast<>(&Type::hm_values))
        .def_property_readonly("hm_indptr", &Type::hm_indptr)
        .def_property_readonly("hm_indices", &Type::hm_indices)
        .def_property_readonly("hm_norm_inf", &Type::hm_norm_inf)
        .def_property_readonly("nb_equations", &Type::nb_equations)
        .def_property_readonly("nb_variables", &Type::nb_variables)
        .def_property_readonly("values", py::overload_cast<>(&Type::values))
        .def_property_readonly("equation_bounds", &Type::equation_bounds)
        .def_property_readonly("variable_bounds", &Type::variable_bounds)
        .def_property_readonly("nb_elements_f", &Type::nb_elements_f)
        .def_property_readonly("nb_elements_g", &Type::nb_elements_g)
        // properties
        .def_property("linear_solver", &Type::linear_solver, &Type::set_linear_solver)
        .def_property("f", &Type::f, &Type::set_f)
        .def_property("nb_threads", &Type::nb_threads, &Type::set_nb_threads)
        .def_property("grainsize", &Type::grainsize, &Type::set_grainsize)
        .def_property("sigma", &Type::sigma, &Type::set_sigma)
        .def_property("hm_diagonal", &Type::hm_diagonal, &Type::set_hm_diagonal)
        .def_property("x", py::overload_cast<>(&Type::x, py::const_), py::overload_cast<Ref<const Vector>>(&Type::set_x, py::const_))
        .def_property("variable_multipliers", py::overload_cast<>(&Type::variable_multipliers, py::const_), py::overload_cast<Ref<const Vector>>(&Type::set_variable_multipliers, py::const_))
        .def_property("equation_multipliers", py::overload_cast<>(&Type::equation_multipliers, py::const_), py::overload_cast<Ref<const Vector>>(&Type::set_equation_multipliers, py::const_))
        // methods
        .def("add_x", py::overload_cast<Ref<const Vector>>(&Type::add_x, py::const_))
        .def("sub_x", py::overload_cast<Ref<const Vector>>(&Type::sub_x, py::const_))
        .def("variable_index", &Type::variable_index, "variable"_a)
        .def("equation_index", &Type::equation_index, "equation"_a)
        .def("clone", &Type::clone)
        .def("remove_inactive_elements", &Type::remove_inactive_elements)
        .def("compute", &Type::compute<true>, "order"_a = 2, py::call_guard<py::gil_scoped_release>())
        .def("hm_add_diagonal", &Type::hm_add_diagonal, "value"_a)
        .def("hm_inv_v", &Type::hm_inv_v, py::call_guard<py::gil_scoped_release>())
        .def("hm_v", &Type::hm_v)
        .def("f_of", [](Type& self, Ref<const Vector> x) {
            self.set_x(x);
            self.compute<false>(0);
            return self.f();
        },
            "x"_a, py::call_guard<py::gil_scoped_release>())
        .def("g_of", [](Type& self, Ref<const Vector> x) {
            self.set_x(x);
            self.compute<false>(0);
            return Vector(self.g());
        },
            "x"_a, py::call_guard<py::gil_scoped_release>())
        .def("df_of", [](Type& self, Ref<const Vector> x) {
            self.set_x(x);
            self.compute<false>(1);
            return Vector(self.df());
        },
            "x"_a, py::call_guard<py::gil_scoped_release>())
        .def("dg_of", [=](Type& self, Ref<const Vector> x) {
            self.set_x(x);
            self.compute<false>(1);
            return csr_matrix(
                std::make_tuple(self.dg_values(), self.dg_indices(), self.dg_indptr()),
                std::make_pair(self.nb_equations(), self.nb_variables()))
                .release();
        },
            "x"_a, py::call_guard<py::gil_scoped_release>())
        .def("hm_of", [=](Type& self, Ref<const Vector> x) {
            self.set_x(x);
            self.compute<false>(2);
            return csr_matrix(
                std::make_tuple(self.hm_values(), self.hm_indices(), self.hm_indptr()),
                std::make_pair(self.nb_variables(), self.nb_variables()))
                .release();
        },
            "x"_a, py::call_guard<py::gil_scoped_release>())
        .def("hm_v_of", [=](Type& self, Ref<const Vector> x, Ref<const Vector> p) {
            self.set_x(x);
            self.compute<false>(2);
            return self.hm_v(p);
        },
            "x"_a, "p"_a, py::call_guard<py::gil_scoped_release>())
        .def("scale", &Type::scale, "factor"_a);
}

} // namespace eqlib
