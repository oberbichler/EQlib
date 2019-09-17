#pragma once

#include "../Define.h"
#include "../System.h"

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
    : m_info_level(0), m_system(system), m_maxiter(100), m_rtol(1e-6)
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

    void compute_f(OptVar* opt, Workspace* wsp, Params* par, Control* cnt)
    {
        m_system->set_x(opt->X);
        m_system->assemble<0>(false);

        opt->F = wsp->ScaleObj * m_system->f();
    }

    void compute_g(OptVar* opt, Workspace* wsp, Params* par, Control* cnt)
    {
    }

    void compute_df(OptVar* opt, Workspace* wsp, Params* par, Control* cnt)
    {
        m_system->set_x(opt->X);
        m_system->assemble<1>(false);

        for (int i = 0; i < m_system->nb_free_dofs(); i++) {
            wsp->DF.val[i] = wsp->ScaleObj * m_system->g(i);
        }
    }

    void compute_dg(OptVar* opt, Workspace* wsp, Params* par, Control* cnt)
    {
    }

    void compute_h(OptVar* opt, Workspace* wsp, Params* par, Control* cnt)
    {
        // Only scale the F part of HM

        m_system->set_x(opt->X);
        m_system->assemble<2>(false);

        Sparse h = m_system->h().transpose();

        int i = 0;
        int j = h.nonZeros() - m_system->nb_free_dofs();

        for (int k = 0; k < h.outerSize(); ++k) {
            for (Sparse::InnerIterator it(h, k); it; ++it){
                if (it.row() == it.col()) {
                    wsp->HM.val[j++] = wsp->ScaleObj * it.value();
                } else {
                    wsp->HM.val[i++] = wsp->ScaleObj * it.value();
                }
            }
        }
    }

    void minimize()
    {
        // setup

        Log::info(1, "==> Minimizing nonlinear system...");
        Log::info(2, "Using WORHP");

        Timer timer;

        OptVar opt;
        Workspace wsp;
        Params par;
        Control cnt;

        CHECK_WORHP_VERSION

        WorhpPreInit(&opt, &wsp, &par, &cnt);

        int status;
        InitParams(&status, &par);

        par.NLPprint = 0;

        ReadParamsNoInit(&status, "worhp.xml", &par);

        if (status == DataError || status == InitError) {
            Log::error("Could not read settings from XML");
            return;;
        }

        opt.n = m_system->nb_free_dofs();
        opt.m = 0;

        wsp.DF.nnz = m_system->nb_free_dofs();
        wsp.DG.nnz = 0;
        wsp.HM.nnz = m_system->h_nb_nonzeros();

        WorhpInit(&opt, &wsp, &par, &cnt);

        if (cnt.status != FirstCall) {
            Log::error("Initialisation failed");
            return;
        }

        for (int i = 0; i < m_system->nb_free_dofs(); i++) {
            opt.X[i] = m_system->dof(i)->delta();
            opt.Lambda[i] = 0.0;
            opt.XL[i] = m_system->dof(i)->lower_bound();
            opt.XU[i] = m_system->dof(i)->upper_bound();
        }

        if (wsp.DF.NeedStructure) {
            for (int i = 0; i < m_system->nb_free_dofs(); i++) {
                wsp.DF.row[i] = i + 1;
            }
        }

        if (wsp.DG.NeedStructure) {
        }

        if (wsp.HM.NeedStructure) {
            int i = 0;
            int j = m_system->h_nb_nonzeros() - m_system->nb_free_dofs();

            Sparse h = m_system->h().transpose();

            for (int k = 0; k < h.outerSize(); ++k) {
                for (Sparse::InnerIterator it(h, k); it; ++it){
                    const int row = it.row() + 1;
                    const int col = it.col() + 1;

                    if (row == col) {
                        wsp.HM.row[j] = row;
                        wsp.HM.col[j] = col;
                        j += 1;
                    } else {
                        wsp.HM.row[i] = row;
                        wsp.HM.col[i] = col;
                        i += 1;
                    }
                }
            }
        }

        while (cnt.status < TerminateSuccess && cnt.status > TerminateError) {
            if (GetUserAction(&cnt, callWorhp)) {
                ::Worhp(&opt, &wsp, &par, &cnt);
            }

            if (GetUserAction(&cnt, iterOutput)) {
                IterationOutput(&opt, &wsp, &par, &cnt);
                DoneUserAction(&cnt, iterOutput);
            }

            if (GetUserAction(&cnt, evalF)) {
                compute_f(&opt, &wsp, &par, &cnt);
                DoneUserAction(&cnt, evalF);
            }

            if (GetUserAction(&cnt, evalG)) {
                compute_g(&opt, &wsp, &par, &cnt);
                DoneUserAction(&cnt, evalG);
            }

            if (GetUserAction(&cnt, evalDF)) {
                compute_df(&opt, &wsp, &par, &cnt);
                DoneUserAction(&cnt, evalDF);
            }

            if (GetUserAction(&cnt, evalDG)) {
                compute_dg(&opt, &wsp, &par, &cnt);
                DoneUserAction(&cnt, evalDG);
            }

            if (GetUserAction(&cnt, evalHM)) {
                compute_h(&opt, &wsp, &par, &cnt);
                DoneUserAction(&cnt, evalHM);
            }

            if (GetUserAction(&cnt, fidif)) {
                WorhpFidif(&opt, &wsp, &par, &cnt);
            }
        }

        StatusMsg(&opt, &wsp, &par, &cnt);

        WorhpFree(&opt, &wsp, &par, &cnt);

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