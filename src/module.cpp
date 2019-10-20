#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <EQlib/Constraint.h>
#include <EQlib/Element.h>
#include <EQlib/Equation.h>
#include <EQlib/Log.h>
#include <EQlib/Node.h>
#include <EQlib/Objective.h>
#include <EQlib/Parameter.h>
#include <EQlib/Problem.h>
#include <EQlib/Point.h>
#include <EQlib/System.h>
#include <EQlib/Variable.h>

#include <EQlib/optimizer/GradientDescent.h>
#ifdef EQLIB_USE_IPOPT
#include <EQlib/optimizer/IPOpt.h>
#endif
#include <EQlib/optimizer/NewtonRaphson.h>
#ifdef EQLIB_USE_WORHP
#include <EQlib/optimizer/Worhp.h>
#endif

#include <EQlib/solver/LBfgsSolver.h>
#include <EQlib/solver/LevenbergMarquardtSolver.h>
#include <EQlib/solver/NewtonDescentSolver.h>
#ifdef EQLIB_USE_WORHP
#include <EQlib/solver/WorhpSolver.h>
#endif

#include <EQlib/Elements/BoundaryConstraint.h>
#include <EQlib/Elements/EqualityConstraint.h>
#include <EQlib/Elements/EqualSubdivisionConstraint.h>
#include <EQlib/Elements/LengthConstraint.h>
#include <EQlib/Elements/NodalEquilibrium.h>

#include <EQlib/Version.h>

PYBIND11_MODULE(EQlib, m) {
    m.doc() = "EQlib by Thomas Oberbichler";
    m.attr("__author__") = "Thomas Oberbichler";
    m.attr("__copyright__") = "Copyright (c) 2018-2019, Thomas Oberbichler";
    m.attr("__version__") = EQlib::version();
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

    // Element
    EQlib::Element::register_python(m);

    // Log
    EQlib::Log::register_python(m);

    // Node
    EQlib::Node::register_python(m);

    // Point
    EQlib::Point::register_python(m);

    // Parameter
    EQlib::Parameter::register_python(m);

    // System
    EQlib::System<false>::register_python(m);

    // SymmetricSystem
    EQlib::System<true>::register_python(m);

    // Timer
    EQlib::Timer::register_python(m);


    // --- solver

    auto solver = m.def_submodule("solver");

    // LBfgsSolver
    EQlib::LBfgsSolver::register_python(solver);

    // LevenbergMarquardtSolver
    EQlib::LevenbergMarquardtSolver::register_python(solver);

    // NewtonDescentSolver
    EQlib::NewtonDescentSolver::register_python(solver);

    // Worhp
    #ifdef EQLIB_USE_WORHP
    EQlib::WorhpSolver::register_python(solver);
    #endif


    // --- optimizer

    auto optimizer = m.def_submodule("optimizer");

    // GradientDescent
    EQlib::GradientDescent::register_python(optimizer);

    // IPOpt
    #ifdef EQLIB_USE_IPOPT
    EQlib::IPOpt::register_python(optimizer);
    #endif

    // NewtonRaphson
    EQlib::NewtonRaphson::register_python(optimizer);

    // Worhp
    #ifdef EQLIB_USE_WORHP
    EQlib::Worhp::register_python(optimizer);
    #endif


    // --- Elements

    // BoundaryConstraint
    EQlib::BoundaryConstraint::register_python(m);

    // EqualityConstraint
    EQlib::EqualityConstraint::register_python(m);

    // EqualSubdivisionConstraint
    EQlib::EqualSubdivisionConstraint::register_python(m);

    // LengthConstraint
    EQlib::LengthConstraint::register_python(m);

    // NodalEquilibrium
    EQlib::NodalEquilibrium::register_python(m);


    // Equation
    EQlib::Equation::register_python(m);

    // Variable
    EQlib::Variable::register_python(m);

    // Objective
    EQlib::Objective::register_python(m);

    // Constraint
    EQlib::Constraint::register_python(m);

    // Problem
    EQlib::Problem::register_python(m);
}
