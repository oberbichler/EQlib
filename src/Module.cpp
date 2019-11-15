#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <eqlib/Constraint.h>
#include <eqlib/Equation.h>
#include <eqlib/LambdaConstraint.h>
#include <eqlib/LambdaObjective.h>
#include <eqlib/Log.h>
#include <eqlib/Objective.h>
#include <eqlib/Problem.h>
#include <eqlib/Node.h>
#include <eqlib/Variable.h>

#include <eqlib/GradientDescent.h>
#include <eqlib/LBfgs.h>
#include <eqlib/LevenbergMarquardt.h>
#include <eqlib/NewtonRaphson.h>
#ifdef EQLIB_USE_WORHP
#include <eqlib/Worhp.h>
#endif

#include <eqlib/objectives/IgaPointSupport.h>
#include <eqlib/objectives/IgaPointSupportAD.h>
#include <eqlib/objectives/IgaShell3P.h>
#include <eqlib/objectives/IgaShell3PAD.h>
#include <eqlib/objectives/IgaShell3PRefAD.h>
#include <eqlib/objectives/IgaShell3PLoadAD.h>
#include <eqlib/objectives/IgaShell3PLoadRefAD.h>

#include <eqlib/Info.h>

PYBIND11_MODULE(eqlib, m) {
    m.doc() = "eqlib by Thomas Oberbichler";
    m.attr("__author__") = "Thomas Oberbichler";
    m.attr("__copyright__") = "Copyright (c) 2018-2019, Thomas Oberbichler";
    m.attr("__version__") = eqlib::Info::version();
    m.attr("__email__") = "thomas.oberbichler@gmail.com";
    m.attr("__status__") = "Development";

    namespace py = pybind11;
    using namespace pybind11::literals;

    m.attr("_GIT_COMMIT_HASH") = eqlib::Info::git_commit_hash();
    m.attr("_USE_BLAS") = eqlib::Info::use_blas();
    m.attr("_USE_MKL") = eqlib::Info::use_mkl();


    // --- core

    // Equation
    eqlib::Equation::register_python(m);

    // Variable
    eqlib::Variable::register_python(m);


    // Constraint
    eqlib::Constraint::register_python(m);

    // Objective
    eqlib::Objective::register_python(m);

    // LambdaConstraint
    eqlib::LambdaConstraint::register_python(m);

    // LambdaObjective
    eqlib::LambdaObjective::register_python(m);


    // Problem
    eqlib::Problem::register_python(m);


    // Log
    eqlib::Log::register_python(m);

    // Node
    eqlib::Node::register_python(m);

    // Timer
    eqlib::Timer::register_python(m);


    // --- solver

    // GradientDescent
    eqlib::GradientDescent::register_python(m);

    // LBfgs
    eqlib::LBfgs::register_python(m);

    // LevenbergMarquardt
    eqlib::LevenbergMarquardt::register_python(m);

    // NewtonRaphson
    eqlib::NewtonRaphson::register_python(m);

    // Worhp
    #ifdef EQLIB_USE_WORHP
    eqlib::Worhp::register_python(m);
    #endif


    // --- objectives

    // IgaPointSupport
    eqlib::IgaPointSupport::register_python(m);

    // IgaShell3P
    eqlib::IgaShell3P::register_python(m);

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