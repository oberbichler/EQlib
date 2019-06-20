#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <EQlib/Dof.h>
#include <EQlib/Element.h>
#include <EQlib/IPOpt.h>
#include <EQlib/LBfgs.h>
#include <EQlib/LevenbergMarquardt.h>
#include <EQlib/Log.h>
#include <EQlib/NewtonDescent.h>
#include <EQlib/Node.h>
#include <EQlib/Parameter.h>
#include <EQlib/System.h>

PYBIND11_MODULE(EQlib, m) {
    m.doc() = "EQlib by Thomas Oberbichler";
    m.attr("__author__") = "Thomas Oberbichler";
    m.attr("__copyright__") = "Copyright (c) 2018-2019, Thomas Oberbichler";
    m.attr("__version__") = EQLIB_VERSION;
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

    // Dof
    EQlib::Dof::register_python(m);

    // Element
    EQlib::Element::register_python(m);

    // Log
    EQlib::Log::register_python(m);

    // Node
    EQlib::Node::register_python(m);

    // Parameter
    EQlib::Parameter::register_python(m);

    // System
    EQlib::System<false>::register_python(m);

    // SymmetricSystem
    EQlib::System<true>::register_python(m);

    // Timer
    EQlib::Timer::register_python(m);

    // --- Solver

    // IPOpt
    EQlib::IPOpt::register_python(m);

    // LBfgs
    EQlib::LBfgs::register_python(m);

    // LevenbergMarquardt
    EQlib::LevenbergMarquardt::register_python(m);

    // NewtonDescent
    EQlib::NewtonDescent::register_python(m);
}