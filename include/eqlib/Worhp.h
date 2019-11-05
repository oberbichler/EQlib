#pragma once

#include "../Define.h"
#include "../Log.h"
#include "../Timer.h"
#include "../Problem.h"

#include <worhp/worhp.h>

#include <cassert>

namespace eqlib {

class Worhp
{
public:     // types

private:    // variables
    Pointer<Problem> m_problem;
    int m_info_level;
    int m_maxiter;
    double m_rtol;

public:     // constructors
    Worhp(Pointer<Problem> problem)
    : m_info_level(0), m_problem(problem), m_maxiter(100), m_rtol(1e-6)
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
        m_problem->set_sigma(wsp->ScaleObj);
        m_problem->set_x(opt->X);
        m_problem->set_equation_multipliers(opt->Mu);

        m_problem->compute();

        opt->F = m_problem->f();
    }

    void compute_g(OptVar* opt, Workspace* wsp, Params* par, Control* cnt)
    {
        m_problem->set_sigma(wsp->ScaleObj);
        m_problem->set_x(opt->X);
        m_problem->set_equation_multipliers(opt->Mu);

        m_problem->compute();

        for (int i = 0; i < m_problem->nb_equations(); i++) {
            opt->G[i] = m_problem->g(i);
        }
    }

    void compute_df(OptVar* opt, Workspace* wsp, Params* par, Control* cnt)
    {
        m_problem->set_sigma(wsp->ScaleObj);
        m_problem->set_x(opt->X);
        m_problem->set_equation_multipliers(opt->Mu);

        m_problem->compute();

        for (int i = 0; i < m_problem->nb_variables(); i++) {
            wsp->DF.val[i] = m_problem->df(i);
        }
    }

    void compute_dg(OptVar* opt, Workspace* wsp, Params* par, Control* cnt)
    {
        m_problem->set_sigma(wsp->ScaleObj);
        m_problem->set_x(opt->X);
        m_problem->set_equation_multipliers(opt->Mu);

        m_problem->compute();

        const Sparse& dg = m_problem->dg();

        int i = 0;

        for (int k = 0; k < dg.outerSize(); ++k) {
            for (Sparse::InnerIterator it(dg, k); it; ++it){
                wsp->DG.val[i++] = it.value();
            }
        }
    }

    void compute_hl(OptVar* opt, Workspace* wsp, Params* par, Control* cnt)
    {
        m_problem->set_sigma(wsp->ScaleObj);
        m_problem->set_x(opt->X);
        m_problem->set_equation_multipliers(opt->Mu);

        m_problem->compute();

        const Sparse hl = m_problem->hl();

        int i = 0;
        int j = hl.nonZeros() - m_problem->nb_variables();

        for (int k = 0; k < hl.outerSize(); ++k) {
            for (Sparse::InnerIterator it(hl, k); it; ++it){
                if (it.row() == it.col()) {
                    wsp->HM.val[j++] = it.value();
                } else {
                    wsp->HM.val[i++] = it.value();
                }
            }
        }
    }

    void run()
    {
        // setup

        Log::info(1, "==> Minimizing nonlinear problem...");
        Log::info(2, "Using WORHP");

        Timer timer;

        OptVar opt;
        Workspace wsp;
        Params par;
        Control cnt;

        // const auto print = [](int mode, const char* message)
        // {
        //     switch (mode) {
        //     case WORHP_PRINT_WARNING:
        //         Log::warn(message);
        //         break;
        //     case WORHP_PRINT_ERROR:
        //         Log::error(message);
        //         break;
        //     default:
        //         Log::info(message);
        //         break;
        //     }
        // };

        // SetWorhpPrint(print);

        CHECK_WORHP_VERSION

        WorhpPreInit(&opt, &wsp, &par, &cnt);

        int status;
        InitParams(&status, &par);

        par.NLPprint = 0;

        Log::info(2, "Read settings from XML...");

        ReadParamsNoInit(&status, "worhp.xml", &par);

        if (status == DataError || status == InitError) {
            Log::error("Could not read settings from XML");
            return;;
        }

        Log::info(2, "Setup optimization problem...");

        opt.n = m_problem->nb_variables();
        opt.m = m_problem->nb_equations();

        const int nnz = (m_problem->hl().nonZeros() + opt.n) / 2;

        wsp.DF.nnz = m_problem->nb_variables();
        wsp.DG.nnz = m_problem->dg().nonZeros();
        wsp.HM.nnz = m_problem->hl().nonZeros();

        Log::info(2, "Initialize Worhp...");

        WorhpInit(&opt, &wsp, &par, &cnt);

        if (cnt.status != FirstCall) {
            Log::error("Initialisation failed");
            return;
        }

        Log::info(2, "Apply variable data...");

        for (int i = 0; i < m_problem->nb_variables(); i++) {
            opt.X[i] = m_problem->variable(i)->value();
            opt.XL[i] = m_problem->variable(i)->lower_bound();
            opt.XU[i] = m_problem->variable(i)->upper_bound();
            opt.Lambda[i] = m_problem->variable(i)->multiplier();
        }

        Log::info(2, "Apply equation data...");

        for (int i = 0; i < m_problem->nb_equations(); i++) {
            opt.GL[i] = m_problem->equation(i)->lower_bound();
            opt.GU[i] = m_problem->equation(i)->upper_bound();
            opt.Mu[i] = m_problem->equation(i)->multiplier();
        }

        if (wsp.DF.NeedStructure) {
            Log::info(2, "Read structure of df...");

            for (int i = 0; i < m_problem->nb_variables(); i++) {
                wsp.DF.row[i] = i + 1;
            }
        }

        if (wsp.DG.NeedStructure) {
            Log::info(2, "Read structure of dg...");

            const Sparse dg = m_problem->dg();

            int i = 0;

            for (int k = 0; k < dg.outerSize(); ++k) {
                for (Sparse::InnerIterator it(dg, k); it; ++it){
                    const int row = it.row() + 1;
                    const int col = it.col() + 1;

                    wsp.DG.row[i] = row;
                    wsp.DG.col[i] = col;

                    i += 1;
                }
            }
        }

        if (wsp.HM.NeedStructure) {
            Log::info(2, "Read structure of hl...");

            int i = 0;
            int j = m_problem->hl().nonZeros() - m_problem->nb_variables();

            const Sparse hl = m_problem->hl();

            for (int k = 0; k < hl.outerSize(); ++k) {
                for (Sparse::InnerIterator it(hl, k); it; ++it){
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
                compute_hl(&opt, &wsp, &par, &cnt);
                DoneUserAction(&cnt, evalHM);
            }

            if (GetUserAction(&cnt, fidif)) {
                WorhpFidif(&opt, &wsp, &par, &cnt);
            }
        }

        StatusMsg(&opt, &wsp, &par, &cnt);

        m_problem->set_sigma(wsp.ScaleObj);
        m_problem->set_x(opt.X);
        m_problem->set_equation_multipliers(opt.Mu);

        WorhpFree(&opt, &wsp, &par, &cnt);

        Log::info(1, "System minimized in {:.3f} sec", timer.ellapsed());
    }

    double rtol() const
    {
        return m_rtol;
    }

    Pointer<Problem> problem() const
    {
        return m_problem;
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
            .def(py::init<Pointer<eqlib::Problem>>(), "problem"_a)
            // read-only properties
            .def_property_readonly("problem", &Type::problem)
            // properties
            .def_property("info_level", &Type::info_level,
                &Type::set_info_level)
            .def_property("maxiter", &Type::maxiter,
                &Type::set_maxiter)
            .def_property("rtol", &Type::rtol, &Type::set_rtol)
            // methods
            .def("run", &Type::run)
        ;
    }
};

} // namespace eqlib