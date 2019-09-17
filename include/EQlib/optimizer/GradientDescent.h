#pragma once

#include "../Define.h"
#include "../Settings.h"
#include "../Timer.h"
#include "../Problem.h"

namespace EQlib {

class GradientDescent
{
private:    // types

private:    // members
    Pointer<Problem> m_problem;
    index m_iterations;
    index m_maxiter;
    index m_fevals;
    index m_gevals;
    double m_rnorm;
    double m_xnorm;
    double m_rtol;
    double m_xtol;

private:    // methods
    double linesearch_armijo(Ref<const Vector> x,
        Ref<const Vector> search_dir, const double alpha_init = 1.0)
    {
        const double c = 0.2;
        const double rho = 0.9;
        double alpha = alpha_init;

        const double f_in = m_problem->f();
        const Vector grad = m_problem->df();

        m_problem->set_x(x + alpha * search_dir);
        m_problem->compute(); //assemble<0>(false);
        m_fevals += 1;
        double f = m_problem->f();

        const double cache = c * grad.dot(search_dir);

        while(f > f_in + alpha * cache) {
            alpha *= rho;

            m_problem->set_x(x + alpha * search_dir);
            m_problem->compute(); //assemble<0>(false);
            m_fevals += 1;
            f = m_problem->f();
        }

        return alpha;
    }

    double linesearch_morethuente(Ref<const Vector> x,
        Ref<const Vector> search_dir, const double alpha_init = 1.0)
    {
        // assume step width
        double ak = alpha_init;

        m_problem->set_x(x);
        m_problem->compute(); //assemble<1>(false);
        m_fevals += 1;
        m_gevals += 1;
        double fval = m_problem->f();
        Vector g = m_problem->df();

        Vector s = search_dir.eval();
        Vector xx = x.eval();

        cvsrch(xx, fval, g, ak, s);

        return ak;
    }

    int cvsrch(Ref<Vector> x, double f, Ref<Vector> g, double &stp,
        Ref<Vector> s)
    {
        // we rewrite this from MIN-LAPACK and some MATLAB code
        int info = 0;
        int infoc = 1;
        const double xtol = 1e-15;
        const double ftol = 1e-4;
        const double gtol = 1e-2;
        const double stpmin = 1e-15;
        const double stpmax = 1e15;
        const double xtrapf = 4;
        const int maxfev = 20;
        int nfev = 0;

        double dginit = g.dot(s);

        if (dginit >= 0.0) {
            // no descent direction
            // TODO: handle this case
            return -1;
        }

        bool brackt      = false;
        bool stage1      = true;

        double finit      = f;
        double dgtest     = ftol * dginit;
        double width      = stpmax - stpmin;
        double width1     = 2 * width;
        Vector wa = x.eval();

        double stx        = 0.0;
        double fx         = finit;
        double dgx        = dginit;
        double sty        = 0.0;
        double fy         = finit;
        double dgy        = dginit;

        double stmin;
        double stmax;

        while (true) {
            // make sure we stay in the interval when setting min/max-step-width
            if (brackt) {
                stmin = std::min(stx, sty);
                stmax = std::max(stx, sty);
            } else {
                stmin = stx;
                stmax = stp + xtrapf * (stp - stx);
            }

            // Force the step to be within the bounds stpmax and stpmin.
            stp = std::max(stp, stpmin);
            stp = std::min(stp, stpmax);

            // Oops, let us return the last reliable values
            if ((brackt && ((stp <= stmin) || (stp >= stmax)))
                || (nfev >= maxfev - 1 ) || (infoc == 0)
                || (brackt && ((stmax - stmin) <= (xtol * stmax)))) {
                stp = stx;
            }

            // test new point
            x = wa + stp * s;

            m_problem->set_x(x);
            m_problem->compute(); //assemble<1>(false);
            m_fevals += 1;
            m_gevals += 1;
            f = m_problem->f();
            g = m_problem->df();
            nfev++;

            double dg = g.dot(s);
            double ftest1 = finit + stp * dgtest;

            // all possible convergence tests
            if ((brackt & ((stp <= stmin) | (stp >= stmax))) | (infoc == 0)) {
                info = 6;
            }

            if ((stp == stpmax) & (f <= ftest1) & (dg <= dgtest)) {
                info = 5;
            }

            if ((stp == stpmin) & ((f > ftest1) | (dg >= dgtest))) {
                info = 4;
            }

            if (nfev >= maxfev) {
                info = 3;
            }

            if (brackt & (stmax - stmin <= xtol * stmax)) {
                info = 2;
            }

            if ((f <= ftest1) & (fabs(dg) <= gtol * (-dginit))) {
                info = 1;
            }

            // terminate when convergence reached
            if (info != 0)
                return -1;

            if (stage1 & (f <= ftest1) & (dg >= std::min(ftol, gtol) * dginit))
                stage1 = false;

            if (stage1 & (f <= fx) & (f > ftest1)) {
                double fm = f - stp * dgtest;
                double fxm = fx - stx * dgtest;
                double fym = fy - sty * dgtest;
                double dgm = dg - dgtest;
                double dgxm = dgx - dgtest;
                double dgym = dgy - dgtest;

                cstep(stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm, brackt,
                    stmin, stmax, infoc);

                fx = fxm + stx * dgtest;
                fy = fym + sty * dgtest;
                dgx = dgxm + dgtest;
                dgy = dgym + dgtest;
            } else {
                // this is ugly and some variables should be moved to the class scope
                cstep(stx, fx, dgx, sty, fy, dgy, stp, f, dg, brackt, stmin,
                    stmax, infoc);
            }

            if (brackt) {
                if (fabs(sty - stx) >= 0.66 * width1) {
                    stp = stx + 0.5 * (sty - stx);
                }
                width1 = width;
                width = fabs(sty - stx);
            }
        }

        return 0;
    }

    int cstep(double& stx, double& fx, double& dx, double& sty, double& fy,
        double& dy, double& stp, double& fp, double& dp, bool& brackt,
        double& stpmin, double& stpmax, int& info)
    {
        info = 0;
        bool bound = false;

        // Check the input parameters for errors.
        if ((brackt & ((stp <= std::min(stx, sty) ) |
            (stp >= std::max(stx, sty)))) |
            (dx * (stp - stx) >= 0.0) | (stpmax < stpmin)) {
            return -1;
        }

        double sgnd = dp * (dx / fabs(dx));

        double stpf = 0;
        double stpc = 0;
        double stpq = 0;

        if (fp > fx) {
            info = 1;
            bound = true;

            double theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp;
            double s = std::max(theta, std::max(dx, dp));
            double gamma = s * sqrt((theta / s) * (theta / s) -
                (dx / s) * (dp / s));

            if (stp < stx) {
                gamma = -gamma;
            }

            double p = (gamma - dx) + theta;
            double q = ((gamma - dx) + gamma) + dp;
            double r = p / q;

            stpc = stx + r * (stp - stx);
            stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2.0) *
                (stp - stx);

            if (fabs(stpc - stx) < fabs(stpq - stx)) {
                stpf = stpc;
            } else {
                stpf = stpc + (stpq - stpc) / 2;
            }

            brackt = true;
        } else if (sgnd < 0.0) {
            info = 2;
            bound = false;

            double theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
            double s = std::max(theta, std::max(dx, dp));
            double gamma = s * sqrt((theta / s) * (theta / s)  - (dx / s) * (dp / s));

            if (stp > stx) {
                gamma = -gamma;
            }

            double p = (gamma - dp) + theta;
            double q = ((gamma - dp) + gamma) + dx;
            double r = p / q;

            stpc = stp + r * (stx - stp);
            stpq = stp + (dp / (dp - dx)) * (stx - stp);

            if (fabs(stpc - stp) > fabs(stpq - stp)) {
                stpf = stpc;
            } else {
                stpf = stpq;
            }

            brackt = true;
        } else if (fabs(dp) < fabs(dx)) {
            info = 3;
            bound = 1;

            double theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
            double s = std::max(theta, std::max( dx, dp));
            double gamma = s * sqrt(std::max(0.0, (theta / s) * (theta / s) -
                (dx / s) * (dp / s)));

            if (stp > stx) {
                gamma = -gamma;
            }

            double p = (gamma - dp) + theta;
            double q = (gamma + (dx - dp)) + gamma;
            double r = p / q;

            if ((r < 0.0) & (gamma != 0.0)) {
                stpc = stp + r * (stx - stp);
            } else if (stp > stx) {
                stpc = stpmax;
            } else {
                stpc = stpmin;
            }

            stpq = stp + (dp / (dp - dx)) * (stx - stp);

            if (brackt) {
                if (fabs(stp - stpc) < fabs(stp - stpq)) {
                    stpf = stpc;
                } else {
                    stpf = stpq;
                }
            } else {
                if (fabs(stp - stpc) > fabs(stp - stpq)) {
                    stpf = stpc;
                } else {
                    stpf = stpq;
                }
            }
        } else {
            info = 4;
            bound = false;

            if (brackt) {
                double theta = 3 * (fp - fy) / (sty - stp) + dy + dp;
                double s = std::max(theta, std::max(dy, dp));
                double gamma = s * sqrt((theta / s) * (theta / s) - (dy / s) * (dp / s));

                if (stp > sty) {
                    gamma = -gamma;
                }

                double p = (gamma - dp) + theta;
                double q = ((gamma - dp) + gamma) + dy;
                double r = p / q;

                stpc = stp + r * (sty - stp);
                stpf = stpc;
            } else if (stp > stx)
                stpf = stpmax;
            else {
                stpf = stpmin;
            }
        }

        if (fp > fx) {
            sty = stp;
            fy = fp;
            dy = dp;
        } else {
            if (sgnd < 0.0) {
                sty = stx;
                fy = fx;
                dy = dx;
            }

            stx = stp;
            fx = fp;
            dx = dp;
        }

        stpf = std::min(stpmax, stpf);
        stpf = std::max(stpmin, stpf);
        stp = stpf;

        if (brackt & bound) {
            if (sty > stx) {
                stp = std::min(stx + 0.66 * (sty - stx), stp);
            } else {
                stp = std::max(stx + 0.66 * (sty - stx), stp);
            }
        }

        return 0;
    }

public:     // constructor
    GradientDescent(Pointer<Problem> problem)
    : m_problem(std::move(problem))
    , m_maxiter(100)
    , m_rtol(1e-6)
    , m_xtol(1e-6)
    {
        if (m_problem->is_constrained()) {
            throw std::runtime_error("Constraints are not supported");
        }
    }

public:     // methods
    index iterations() const noexcept
    {
        return m_iterations;
    }

    index fevals() const noexcept
    {
        return m_fevals;
    }

    index gevals() const noexcept
    {
        return m_gevals;
    }

    index maxiter() const noexcept
    {
        return m_maxiter;
    }

    void set_maxiter(const index value) noexcept
    {
        m_maxiter = value;
    }

    double rnorm() const noexcept
    {
        return m_rnorm;
    }

    double rtol() const noexcept
    {
        return m_rtol;
    }

    void set_rtol(const double value) noexcept
    {
        m_rtol = value;
    }

    double xtol() const noexcept
    {
        return m_xtol;
    }

    void set_xtol(const double value) noexcept
    {
        m_xtol = value;
    }

    void run()
    {
        // setup

        Log::info(1, "==> Minimizing nonlinear problem...");
        Log::info(2, "Using Gradient Descent minimizer");

        Timer timer;

        const auto n = m_problem->nb_variables();

        Vector x = m_problem->x();
        Vector delta(n);
        Vector direction(n);

        m_iterations = 0;
        m_fevals = 0;
        m_gevals = 0;

        while (true) {
            m_problem->compute(); //assemble<1>(false);
            m_fevals += 1;
            m_gevals += 1;

            direction = -m_problem->df();

            m_rnorm = direction.norm();
            
            Log::info(2, "The norm of the residual is {}", m_rnorm);

            if (m_rnorm < m_rtol) {
                break;
            }

            const double rate = linesearch_morethuente(x, direction, 1.0);

            delta = rate * direction;
            
            Log::info(2, "The line search rate is {}", rate);

            m_xnorm = delta.norm();
            
            Log::info(2, "The norm of the step is {}", m_xnorm);

            if (m_xnorm < m_xtol) {
                break;
            }

            x += delta;

            m_problem->set_x(x);

            m_iterations += 1;
            
            if (m_iterations >= m_maxiter) {
                break;
            }
        }

        Log::info(2, "{} iterations", m_iterations);

        Log::info(1, "Problem minimized in {:.3f} sec", timer.ellapsed());
    }

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Type = EQlib::GradientDescent;

        py::class_<Type>(m, "GradientDescent")
            .def(py::init<Pointer<EQlib::Problem>>(), "problem"_a)
            .def("run", &Type::run)
            .def_property("rtol", &Type::rtol, &Type::set_rtol)
            .def_property("maxiter", &Type::maxiter, &Type::set_maxiter)
            .def_property_readonly("iterations", &Type::iterations)
            .def_property_readonly("rnorm", &Type::rnorm)
            .def_property_readonly("fevals", &Type::fevals)
            .def_property_readonly("gevals", &Type::gevals)
        ;
    }
};

} // namespace EQlib