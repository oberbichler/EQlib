#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <EQlib/Dof.h>
#include <EQlib/Log.h>
#include <EQlib/Node.h>
#include <EQlib/Parameter.h>
#include <EQlib/PyElement.h>
#include <EQlib/System.h>

#include <IGAlib/LocationConstraint.h>
#include <IGAlib/Shell3D3P.h>

PYBIND11_MODULE(IGAlib, m) {
    m.doc() = "IGAlib by Thomas Oberbichler";
    m.attr("__author__") = "Thomas Oberbichler";
    m.attr("__copyright__") = "Copyright (c) 2019, Thomas Oberbichler";
    m.attr("__version__") = EQlib_VERSION;
    m.attr("__email__") = "thomas.oberbichler@gmail.com";
    m.attr("__status__") = "Development";

    namespace py = pybind11;
    using namespace pybind11::literals;

#if defined(EIGEN_USE_BLAS)
    m.attr("USE_BLAS") = true;
#else
    m.attr("USE_BLAS") = false;
#endif // EIGEN_USE_BLAS

#if defined(EIGEN_USE_MKL_ALL)
    m.attr("USE_MKL") = true;
#else
    m.attr("USE_MKL") = false;
#endif // EIGEN_USE_MKL_ALL

    // LocationConstraint
    {
        using Type = EQlib::LocationConstraint;
        using Base = EQlib::Element;
        using Holder = std::shared_ptr<Type>;

        py::class_<Type, Base, Holder>(m, "LocationConstraint")
            .def(py::init<std::vector<std::shared_ptr<EQlib::Node>>,
                EQlib::Matrix, EQlib::Vector3D, double>(), "nodes"_a,
                "shape_functions"_a, "target"_a, "penalty"_a=1)
        ;
    }

    // Shell3D3P
    {
        using Type = EQlib::Shell3D3P;
        using Base = EQlib::Element;
        using Holder = std::shared_ptr<Type>;

        py::class_<Type, Base, Holder>(m, "Shell3D3P")
            .def(py::init<std::vector<std::shared_ptr<EQlib::Node>>,
                EQlib::Matrix, double, double, double, double>(), "nodes"_a,
                "shape_functions"_a, "thickness"_a, "young_modulus"_a,
                "poisson_ratio"_a, "weight"_a)
        ;
    }
}