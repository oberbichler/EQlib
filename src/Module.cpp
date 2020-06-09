#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <eqlib/Armijo.h>
#include <eqlib/Constraint.h>
#include <eqlib/Equation.h>
#include <eqlib/LambdaConstraint.h>
#include <eqlib/LambdaObjective.h>
#include <eqlib/Log.h>
#include <eqlib/NewtonRaphson.h>
#include <eqlib/Node.h>
#include <eqlib/Objective.h>
#include <eqlib/Parameter.h>
#include <eqlib/Problem.h>
#include <eqlib/SteepestDecent.h>
#include <eqlib/SparseLU.h>
#include <eqlib/SparseStructure.h>
#include <eqlib/Variable.h>

#include <eqlib/Info.h>

#include <eqlib/objectives/IgaNormalDistanceAD.h>
#include <eqlib/objectives/IgaPointDistance.h>
#include <eqlib/objectives/IgaPointDistanceAD.h>
#include <eqlib/objectives/IgaPointLocation.h>
#include <eqlib/objectives/IgaRotationCouplingAD.h>
#include <eqlib/objectives/IgaShell3PAD.h>

PYBIND11_MODULE(eqlib, m)
{
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.doc() = "eqlib by Thomas Oberbichler";
    m.attr("__author__") = "Thomas Oberbichler";
    m.attr("__copyright__") = "Copyright (c) 2018-2019, Thomas Oberbichler";
    m.attr("__version__") = eqlib::Info::version();
    m.attr("__email__") = "thomas.oberbichler@gmail.com";
    m.attr("__status__") = "Development";

    m.attr("_GIT_COMMIT_HASH") = eqlib::Info::git_commit_hash();
    m.attr("_USE_BLAS") = eqlib::Info::use_blas();
    m.attr("_USE_MKL") = eqlib::Info::use_mkl();

    // --- core

    // Armijo
    eqlib::Armijo::register_python(m);

    // Equation
    eqlib::Equation::register_python(m);

    // Variable
    eqlib::Variable::register_python(m);

    // Parameter
    eqlib::Parameter::register_python(m);

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

    // NewtonRaphson
    eqlib::NewtonRaphson::register_python(m);

    // SteepestDecent
    eqlib::SteepestDecent::register_python(m);

    // LinearSolver
    eqlib::LinearSolver::register_python(m);

    // SimplicialLDLT
    eqlib::SimplicialLDLT::register_python(m);

    // SparseLU
    eqlib::SparseLU::register_python(m);

    #ifdef EQLIB_USE_MKL
    // PardisoLDLT
    eqlib::PardisoLDLT::register_python(m);
    #endif

    // SparseStructure
    eqlib::SparseStructure<double, int, true, true>::register_python(m, "RowMajorSparseStructure");
    eqlib::SparseStructure<double, int, false, true>::register_python(m, "ColMajorSparseStructure");

    // objectives: IgaNormalDistance
    eqlib::IgaNormalDistanceAD::register_python(m);

    // objectives: IgaPointDistance
    eqlib::IgaPointDistance::register_python(m);

    // objectives: IgaPointDistanceAD
    eqlib::IgaPointDistanceAD::register_python(m);

    // objectives: IgaPointLocation
    eqlib::IgaPointLocation::register_python(m);

    // objectives: IgaRotationCouplingAD
    eqlib::IgaRotationCouplingAD::register_python(m);

    // objectives: IgaShell3PAD
    eqlib::IgaShell3PAD::register_python(m);
}