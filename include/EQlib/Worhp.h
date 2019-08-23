#pragma once

#include "Define.h"
#include "System.h"

#include <worhp/worhp.h>

#include <cassert>

namespace EQlib {

class Worhp
{
public:     // types

private:    // variables
    Pointer<System<true>> m_system;
    int m_info_level;
    int m_maxiter;
    double m_rtol;

public:     // constructors
    Worhp(Pointer<System<true>> system)
    : m_info_level(0), m_system(std::move(system)), m_maxiter(100), m_rtol(1e-6)
    { }

public:     // methods
    int info_level() const
    {
        return m_info_level;
    }

    int maxiter() const
    {
        return m_maxiter;
    }

    void minimize()
    {
        // setup

        Log::info(1, "==> Minimizing nonlinear system...");
        Log::info(2, "Using WORHP");

        Timer timer;

        OptVar    opt;
        Workspace wsp;
        Params    par;
        Control   cnt;
    
        // Check Version of library and header files
        CHECK_WORHP_VERSION
    
        // Properly zeros everything, or else the following routines could get confused
        WorhpPreInit(&opt, &wsp, &par, &cnt);
    
        // Uncomment this to get more info on data structures
        //WorhpDiag(&opt, &wsp, &par, &cnt);
    
        /*
        * Parameter initialisation routine that must be called
        * when using ReadParamsNoInit instead of ReadParams.
        */
        int status;
        InitParams(&status, &par);
    
        /*
        * We can now set parameters that may be overruled by those in the
        * parameter file. This is useful for setting a non-default standard
        * parameter value that may still be overwritten.
        */
        par.NLPprint = 1;  // Let's prefer the slim output format
                        // unless the parameter file says differently
    
        /*
        * Parameter XML import routine that does not reset
        * all parameters to default values (InitParams does this)
        */
        ReadParamsNoInit(&status, "worhp.xml", &par);
        if (status == DataError || status == InitError)
        {
            return;// EXIT_FAILURE;
        }
    
        /*
        * WORHP data structure initialisation routine.
        * Calling this routine prior to WORHP is mandatory.
        * Before calling WorhpInit, set the problem and matrix dimensions as
        *
        * opt.n      = number of variables,
        * opt.m      = number of constraints (lin + nonlin, excluding box con's),
        * wsp.DF.nnz = nonzero entries of the objective function gradient,
        * wsp.DG.nnz = nonzero entries of the constraint Jacobian,
        * wsp.HM.nnz = nonzero entries of the Lagrange Hessian.
        *
        * Set nnz to 'WorhpMatrix_Init_Dense' to have WorhpInit allocate and
        * create a dense matrix structure appropriate for the matrix kind and
        * its dimensions. Setting it to its dense dimension achieves the same.
        */
        opt.n = 4;  // This problem has 4 variables
        opt.m = 3;  // and 3 constraints (excluding box constraints)
    
        // All derivatives for this problem have a sparse structure, so
        // set the amount of nonzeros here
        wsp.DF.nnz = 3;
        wsp.DG.nnz = 6;
        wsp.HM.nnz = 1 + opt.n;  // 1 entry on strict lower triangle
                                // plus full diagonal
    
        WorhpInit(&opt, &wsp, &par, &cnt);
        if (cnt.status != FirstCall)
        {
            std::cout << "Main: Initialisation failed." << std::endl;
            return;// EXIT_FAILURE;
        }

        Log::info(1, "System minimized in {:.3f} sec", timer.ellapsed());
    }

    double rtol() const
    {
        return m_rtol;
    }

    Pointer<System<true>> system() const
    {
        return m_system;
    }

    void set_info_level(int value)
    {
        m_info_level = value;
    }

    void set_maxiter(int value)
    {
        m_maxiter = value;
    }

    void set_rtol(double value)
    {
        m_rtol = value;
    }

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Type = Worhp;

        py::class_<Type>(m, "Worhp")
            // constructors
            .def(py::init<Pointer<EQlib::System<true>>>(), "system"_a)
            // read-only properties
            .def_property_readonly("system", &Type::system)
            // properties
            .def_property("info_level", &Type::info_level,
                &Type::set_info_level)
            .def_property("maxiter", &Type::maxiter,
                &Type::set_maxiter)
            .def_property("rtol", &Type::rtol, &Type::set_rtol)
            // methods
            .def("minimize", &Type::minimize)
        ;
    }
};

} // namespace EQlib