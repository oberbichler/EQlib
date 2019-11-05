#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <eqlib/Constraint.h>
#include <eqlib/Equation.h>
#include <eqlib/Log.h>
#include <eqlib/Objective.h>
#include <eqlib/Problem.h>
#include <eqlib/Node.h>
#include <eqlib/Variable.h>

#include <eqlib/GradientDescent.h>
#include <eqlib/LBfgs.h>
#include <eqlib/NewtonRaphson.h>
#ifdef EQLIB_USE_WORHP
#include <eqlib/Worhp.h>
#endif

#include <eqlib/objectives/IgaPointSupportAD.h>
#include <eqlib/objectives/IgaShell3PAD.h>
#include <eqlib/objectives/IgaShell3PRefAD.h>
#include <eqlib/objectives/IgaShell3PLoadAD.h>
#include <eqlib/objectives/IgaShell3PLoadRefAD.h>

#include <eqlib/Version.h>

PYBIND11_MODULE(eqlib, m) {
    m.doc() = "eqlib by Thomas Oberbichler";
    m.attr("__author__") = "Thomas Oberbichler";
    m.attr("__copyright__") = "Copyright (c) 2018-2019, Thomas Oberbichler";
    m.attr("__version__") = eqlib::version();
    m.attr("__email__") = "thomas.oberbichler@gmail.com";
    m.attr("__status__") = "Development";

    namespace py = pybind11;
    using namespace pybind11::literals;

#if defined(GIT_COMMIT_HASH)
    m.attr("GIT_COMMIT_HASH") = GIT_COMMIT_HASH;
#endif // GIT_COMMIT_HASH

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


    // --- core

    // Constraint
    eqlib::Constraint::register_python(m);

    // Equation
    eqlib::Equation::register_python(m);

    // Log
    eqlib::Log::register_python(m);

    // Node
    eqlib::Node::register_python(m);

    // Objective
    eqlib::Objective::register_python(m);

    // Problem
    eqlib::Problem::register_python(m);

    // Timer
    eqlib::Timer::register_python(m);

    // Variable
    eqlib::Variable::register_python(m);


    // --- solver

    // LBfgs
    eqlib::LBfgs::register_python(m);

    // GradientDescent
    eqlib::GradientDescent::register_python(m);

    // NewtonRaphson
    eqlib::NewtonRaphson::register_python(m);

    // Worhp
    #ifdef EQLIB_USE_WORHP
    eqlib::Worhp::register_python(m);
    #endif


    // --- objectives

    // IgaShell3PAD
    eqlib::IgaShell3PAD::register_python(m);

    // IgaShell3PRefAD
    eqlib::IgaShell3PRefAD::register_python(m);

    // IgaShell3PLoadAD
    eqlib::IgaShell3PLoadAD::register_python(m);

    // IgaShell3PLoadRefAD
    eqlib::IgaShell3PLoadRefAD::register_python(m);

    // IgaPointSupportAD
    eqlib::IgaPointSupportAD::register_python(m);
}