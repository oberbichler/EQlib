#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/eval.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

#include <eqlib/constraint.h>
#include <eqlib/objective.h>

namespace py = pybind11;
using namespace py::literals;

template <typename T>
auto bind(py::module &m, const std::string &name)
{
    using namespace eqlib;

    py::class_<T, Pointer<T>> cls(m, name.c_str());

    return cls;
}

template <typename T, typename Trampoline>
auto bind(py::module &m, const std::string &name)
{
    using namespace eqlib;

    py::class_<T, Trampoline, Pointer<T>> cls(m, name.c_str());

    return cls;
}

template <typename T>
auto bind_objective(py::module &m, const std::string &name)
{
    using namespace eqlib;

    py::class_<T, Objective, Pointer<T>> cls(m, name.c_str());

    return cls;
}

template <typename T>
auto bind_constraint(py::module &m, const std::string &name)
{
    using namespace eqlib;

    py::class_<T, Constraint, Pointer<T>> cls(m, name.c_str());

    return cls;
}