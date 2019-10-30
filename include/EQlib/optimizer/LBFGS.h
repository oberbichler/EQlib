#pragma once

#include "../Define.h"
#include "../Log.h"
#include "../Problem.h"
#include "../Timer.h"

#include <Eigen/LU>

namespace EQlib {

class LBFGS
{
private:    // types
    using Type = EQlib::LBFGS;

private:    // members
    Pointer<Problem> m_problem;
    index m_maxiter;
    double m_rnorm;
    double m_rtol;
    double m_xtol;

protected:
    double linesearch_morethuente(Ref<const Vector> x,
        Ref<const Vector> search_dir, const double alpha_init = 1.0)
    {
        // assume step width
        double ak = alpha_init;

        m_problem->set_x(x);
        m_problem->compute(1);
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

        bool bracket      = false;
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
            if (bracket) {
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
            if ((bracket && ((stp <= stmin) || (stp >= stmax)))
                || (nfev >= maxfev - 1 ) || (infoc == 0)
                || (bracket && ((stmax - stmin) <= (xtol * stmax)))) {
                stp = stx;
            }

            // test new point
            x = wa + stp * s;

            m_problem->set_x(x);
            m_problem->compute(1);
            f = m_problem->f();
            g = m_problem->df();
            nfev++;

            double dg = g.dot(s);
            double ftest1 = finit + stp * dgtest;

            // all possible convergence tests
            if ((bracket & ((stp <= stmin) | (stp >= stmax))) | (infoc == 0)) {
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

            if (bracket & (stmax - stmin <= xtol * stmax)) {
                info = 2;
            }

            if ((f <= ftest1) & (fabs(dg) <= gtol * (-dginit))) {
                info = 1;
            }

            // terminate when convergence reached
            if (info != 0) {
                return -1;
            }

            if (stage1 & (f <= ftest1) & (dg >= std::min(ftol, gtol) * dginit))
                stage1 = false;

            if (stage1 & (f <= fx) & (f > ftest1)) {
                double fm = f - stp * dgtest;
                double fxm = fx - stx * dgtest;
                double fym = fy - sty * dgtest;
                double dgm = dg - dgtest;
                double dgxm = dgx - dgtest;
                double dgym = dgy - dgtest;

                cstep(stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm, bracket,
                    stmin, stmax, infoc);

                fx = fxm + stx * dgtest;
                fy = fym + sty * dgtest;
                dgx = dgxm + dgtest;
                dgy = dgym + dgtest;
            } else {
                // this is ugly and some variables should be moved to the class scope
                cstep(stx, fx, dgx, sty, fy, dgy, stp, f, dg, bracket, stmin,
                    stmax, infoc);
            }

            if (bracket) {
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
        double& dy, double& stp, double& fp, double& dp, bool& bracket,
        double& stpmin, double& stpmax, int& info)
    {
        info = 0;
        bool bound = false;

        // Check the input parameters for errors.
        if ((bracket & ((stp <= std::min(stx, sty) ) |
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

            bracket = true;
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

            bracket = true;
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

            if (bracket) {
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

            if (bracket) {
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

        if (bracket & bound) {
            if (sty > stx) {
                stp = std::min(stx + 0.66 * (sty - stx), stp);
            } else {
                stp = std::max(stx + 0.66 * (sty - stx), stp);
            }
        }

        return 0;
    }

public:     // constructor
    LBFGS(Pointer<Problem> system) : m_problem(system), m_maxiter(100),
        m_rtol(1e-6), m_xtol(1e-6)
    {
    }

public:     // methods
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
        Timer timer;

        const size_t m = 10;
        const size_t nb_variables = m_problem->nb_variables();

        Matrix sVector = Matrix::Zero(nb_variables, m);
        Matrix yVector = Matrix::Zero(nb_variables, m);
        Vector alpha = Vector::Zero(m);
        Vector grad(nb_variables);
        Vector q(nb_variables);
        Vector grad_old(nb_variables);
        Vector s(nb_variables);
        Vector y(nb_variables);;

        m_problem->compute(1);
        Log::info(2, "The current value is {}", m_problem->f());

        Vector x0 = m_problem->x();

        grad = m_problem->df();

        Vector x_old = x0;

        int iter = 0;
        double H0k = 1;

        for (int iteration = 0; iteration < m_maxiter; iteration++) {
            const double relativeEpsilon = 0.0001 * std::max(1.0, x0.norm());

            if (grad.norm() < m_rtol) {
                break;
            }

            //Algorithm 7.4 (L-BFGS two-loop recursion)
            q = grad;

            std::cout << "grad: " << grad.norm() << std::endl;

            const int k = std::min<int>(m, iter);

            // for i = k − 1, k − 2, . . . , k − m§
            for (int i = k - 1; i >= 0; i--) {
                Vector s_i = sVector.col(i);
                Vector y_i = yVector.col(i);
                // alpha_i <- rho_i*s_i^T*q
                const double rho = 1.0 / s_i.dot(y_i);
                alpha(i) = rho * s_i.dot(q);
                // q <- q - alpha_i*y_i
                q = q - alpha(i) * yVector.col(i);
            }
            // r <- H_k^0*q
            q = H0k * q;
            //for i k − m, k − m + 1, . . . , k − 1
            for (int i = 0; i < k; i++) {
                // beta <- rho_i * y_i^T * r
                Vector s_i = sVector.col(i);
                Vector y_i = yVector.col(i);
                const double rho = 1.0 / s_i.dot(y_i);
                const double beta = rho * y_i.dot(q);
                // r <- r + s_i * ( alpha_i - beta)
                q = q + sVector.col(i) * (alpha(i) - beta);
            }
            // stop with result "H_k*f_f'=q"

            // any issues with the descent direction ?
            double descent = -grad.dot(q);
            double alpha_init =  1.0 / grad.norm();

            if (descent > -0.0001 * relativeEpsilon) {
                q = -1 * grad;
                iter = 0;
                alpha_init = 1.0;
            }

            grad_old = m_problem->df();

            // find steplength
            const double rate = linesearch_morethuente(x0, -q, alpha_init);

            // update guess
            x0 = x0 - rate * q;

            m_problem->set_x(x0);
            m_problem->compute(1);

            grad = m_problem->df();

            Log::info(2, "The current value is {}", m_problem->f());
            Log::info(2, "The norm of the residual is {}", grad.norm());

            s = x0 - x_old;
            y = grad - grad_old;

            // update the history
            if (iter < m) {
                sVector.col(iter) = s;
                yVector.col(iter) = y;
            } else {
                sVector.leftCols(m - 1) = sVector.rightCols(m - 1).eval();
                sVector.rightCols(1) = s;
                yVector.leftCols(m - 1) = yVector.rightCols(m - 1).eval();
                yVector.rightCols(1) = y;
            }

            // update the scaling factor
            H0k = y.dot(s) / y.dot(y);

            x_old = x0;

            iter++;
        }

        Log::info(1, "System solved in {:.3f} sec", timer.ellapsed());
    }

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        py::class_<Type>(m, "LBFGS")
            // constructors
            .def(py::init<Pointer<EQlib::Problem>>(), "problem"_a)
            // methods
            .def("run", &Type::run)
            // properties
            .def_property("maxiter", &Type::maxiter, &Type::set_maxiter)
            .def_property("rtol", &Type::rtol, &Type::set_rtol)
            // read-only properties
            // .def_property_readonly("iterations", &Type::iterations)
            // .def_property_readonly("rnorm", &Type::rnorm)
            // .def_property_readonly("fevals", &Type::fevals)
            // .def_property_readonly("gevals", &Type::gevals)
            // .def_property_readonly("hevals", &Type::hevals)
        ;
    }
};

} // namespace EQlib