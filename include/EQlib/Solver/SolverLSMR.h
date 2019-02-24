#pragma once

#include "Define.h"
#include "Solver.h"

#include <iosfwd>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

namespace EQlib {

class LsmrSolverBase : public Solver
{
    static inline double
    sqr(
        const double x) noexcept
    {
        return x * x;
    }

    static inline void
    daxpy(
        size_t n,
        double alpha,
        const double* x,
        double* y)
    {
        const double* xend = x+n;
        while ( x!=xend ) {
            *y++ += alpha * *x++;
        }
    }

public:
    LsmrSolverBase()
    : m_eps(1e-16)
    , m_atol(1e-6)
    , m_btol(1e-6)
    , m_conlim(1.0 / (10 * std::sqrt(m_eps)))
    , m_itnlim(10)
    , m_nout(nullptr)
    , m_istop(0)
    , m_itn(0)
    , m_norm_a(0.0)
    , m_cond_a(0.0)
    , m_norm_r(0.0)
    , m_norm_ar(0.0)
    , m_norm_x(0.0)
    , m_norm_b(0.0)
    , m_dxmax(0.0)
    , m_maxdx(0)
    , m_damp(0.0)
    , m_damped(false)
    , m_local_size(0)
    { }

    virtual ~LsmrSolverBase() { }

    virtual void a_prod_1(Ref<const Vector> x, Ref<Vector> y) const = 0;

    virtual void a_prod_2(Ref<Vector> x, Ref<const Vector> y) const = 0;

    double norm(const double a, const double b) const
    {
        const double scale = std::abs(a) + std::abs(b);
        const double zero = 0.0;

        if (scale == zero) {
            return zero;
        }

        const double sa = a / scale;
        const double sb = b / scale;

        return scale * std::sqrt(sa * sa + sb * sb);
    }

    double norm(size_t n, const double* x) const
    {
        double magnitudeOfLargestElement = 0.0;

        double sumOfSquaresScaled = 1.0;

        for (size_t i = 0; i < n; i++) {
            if (x[i] != 0.0) {
                double dx = x[i];
                const double absxi = std::abs(dx);

                if ( magnitudeOfLargestElement < absxi ) {
                    // rescale the sum to the range of the new element
                    dx = magnitudeOfLargestElement / absxi;
                    sumOfSquaresScaled = sumOfSquaresScaled * (dx * dx) + 1.0;
                    magnitudeOfLargestElement = absxi;
                } else {
                    // rescale the new element to the range of the sum
                    dx = absxi / magnitudeOfLargestElement;
                    sumOfSquaresScaled += dx * dx;
                }
            }
        }

        const double norm = magnitudeOfLargestElement * std::sqrt(sumOfSquaresScaled);

        return norm;
    }

    void scale(size_t n, double factor, double *x) const
    {
        double* xend = x + n;
        while (x != xend) {
            *x++ *= factor;
        }
    }

    void set_local_size(size_t n)
    {
        m_local_size = n;
    }

    void set_tolerance_a(const double value)
    {
        m_atol = value;
    }

    void set_tolerance_b(const double value)
    {
        m_btol = value;
    }

    void set_upper_limit_on_conditional(const double value)
    {
        m_conlim = value;
    }

    void set_epsilon(const double value)
    {
        m_eps = value;
    }

    void set_damp(const double value)
    {
        m_damp = value;
    }

    void set_maximum_number_of_iterations(size_t value)
    {
        m_itnlim = value;
    }

    void set_output_stream(std::ostream & os)
    {
        m_nout = &os;
    }

    size_t stopping_reason() const
    {
        return m_istop;
    }

    std::string stopping_reason_message() const
    {
        switch (m_istop) {
        case 0:
            return "The exact solution is  x = 0";
        case 1:
            return "Ax - b is small enough, given m_atol, m_btol";
        case 2:
            return "The least-squares solution is good enough, given m_atol";
        case 3:
            return "The estimate of cond(Abar) has exceeded m_conlim";
        case 4:
            return "Ax - b is small enough for this machine";
        case 5:
            return "The LS solution is good enough for this machine";
        case 6:
            return "Cond(Abar) seems to be too large for this machine";
        case 7:
            return "The iteration limit has been reached";
        default:
            return "Error. Unknown stopping reason";
        }
    }

    size_t number_of_iterations_performed() const
    {
        return m_itn;
    }

    double frobenius_norm_estimate_of_abar() const
    {
        return m_norm_a;
    }

    double condition_number_estimate_of_abar() const
    {
        return m_cond_a;
    }

    double final_estimate_of_norm_rbar() const
    {
        return m_norm_r;
    }

    double final_estimate_of_norm_of_residuals() const
    {
        return m_norm_ar;
    }

    double final_estimate_of_norm_of_x() const
    {
        return m_norm_x;
    }

    void solve(Ref<const Vector> vb, Ref<Vector> vx) override
    {
        const size_t m = static_cast<size_t>(vb.size());
        const size_t n = static_cast<size_t>(vx.size());
        const double* b = vb.data();
        double* x = vx.data();

        const double zero = 0.0;
        const double one = 1.0;

        double test1;
        double test2;

        // Initialize.

        size_t local_vecs = std::min(m_local_size, std::min(m, n));

        if (m_nout) {
            (*m_nout) << " Enter LSMR. Least-squares solution of  Ax = b\n" << std::endl;
            (*m_nout) << " The matrix  A  has " << m << " rows   and " << n << " columns" << std::endl;
            (*m_nout) << " damp = " << m_damp << std::endl;
            (*m_nout) << " atol = " << m_atol << ", conlim = " << m_conlim << std::endl;
            (*m_nout) << " btol = " << m_btol << ", itnlim = " << m_itnlim << std::endl;
            (*m_nout) << " local_size (no. of vectors for local reorthogonalization) = " << m_local_size << std::endl;
        }

        int pfreq = 20;
        int pcount = 0;
        m_damped = (m_damp > zero);

        std::vector<double> workBuffer(m + 5 * n + n * local_vecs);
        double* u = &workBuffer[0];
        double* v = u + m;
        double* w = v + n;
        double* h = w + n;
        double* hbar = h + n;
        double* local_v = hbar + n;

        //-------------------------------------------------------------------
        //  Set up the first vectors u and v for the bidiagonalization.
        //  These satisfy  beta*u = b,  alpha*v = A(transpose)*u.
        //-------------------------------------------------------------------
        std::copy(b, b + m, u);
        std::fill(v, v + n, zero);
        std::fill(w, v + n, zero);
        std::fill(x, x + n, zero);

        scale(m, -1.0, u);
        a_prod_1(Eigen::Map<Vector>(x, n), Eigen::Map<Vector>(u, m));
        scale(m, -1.0, u);

        double alpha = zero;

        double beta = norm(m, u);

        if (beta > zero) {
            scale(m, one / beta, u);
            a_prod_2(Eigen::Map<Vector>(v, n), Eigen::Map<Vector>(u, m)); // v = A'*u
            alpha = norm(n, v);
        }

        if (alpha > zero) {
            scale(n, one / alpha, v);
            std::copy(v, v + n, w);
        }

        m_norm_ar = alpha * beta;

        if (m_norm_ar == zero) {
            termination_print_out();
            return;
        }

        // Initialization for local reorthogonalization.
        bool local_ortho = false;
        bool local_v_queue_full = false;
        size_t local_pointer = 0;
        if (local_vecs > 0) {
            local_ortho = true;
            std::copy(v, v + n, local_v);
        }

        // Initialize variables for 1st iteration.
        m_itn = 0;
        double zetabar = alpha * beta;
        double alphabar = alpha;
        double rho = one;
        double rhobar = one;
        double cbar = one;
        double sbar = zero;

        std::copy(v, v + n, h);
        std::fill(hbar, hbar + n, zero);

        // Initialize variables for estimation of ||r||.

        double betadd = beta;
        double betad = zero;
        double rhodold = one;
        double tautildeold = zero;
        double thetatilde = zero;
        double zeta = zero;
        double d = zero;

        // Initialize variables for estimation of ||A|| and cond(A).

        double norm_a2 = alpha * alpha;
        double maxrbar = zero;
        double minrbar = 1e+100;

        // Items for use in stopping rules.
        m_norm_b = beta;
        m_istop = 0;
        double ctol = zero;

        if (m_conlim > zero) {
            ctol = one / m_conlim;
        }
        m_norm_r  = beta;

        if (m_nout) {
            if (m_damped) {
                (*m_nout) << "   m_Itn       x(1)           norm rbar    Abar'rbar"
                    " Compatible    LS    norm Abar cond Abar\n";
            } else {
                (*m_nout) << "   m_Itn       x(1)            norm r         A'r   "
                    " Compatible    LS      norm A    cond A\n";
            }

            test1 = one;
            test2 = alpha / beta;

            (*m_nout) << m_itn << ", " << x[0] << ", " << m_norm_r << ", " << m_norm_a << ", " << test1 << ", " << test2 << std::endl;
        }

        // Main iteration loop
        do {
            m_itn++;

            //----------------------------------------------------------------
            //  Perform the next step of the bidiagonalization to obtain the
            //  next beta, u, alpha, v.  These satisfy
            //      beta*u = A*v  - alpha*u,
            //     alpha*v = A'*u -  beta*v.
            //----------------------------------------------------------------
            scale(m, -alpha, u);

            a_prod_1(Eigen::Map<Vector>(v, n), Eigen::Map<Vector>(u, m)); // u = u + A * v

            beta = norm(m, u);

            if (beta > zero) {
                scale(m, one / beta, u);
                if (local_ortho) {
                    if (local_pointer + 1 < local_vecs) {
                        local_pointer = local_pointer + 1;
                    } else {
                        local_pointer = 0;
                        local_v_queue_full = true;
                    }
                    std::copy(v, v + n, local_v+local_pointer * n);
                }

                scale(n, -beta, v);
                a_prod_2(Eigen::Map<Vector>(v, n), Eigen::Map<Vector>(u, m)); // v = A'*u

                if (local_ortho) {
                    size_t localOrthoLimit = local_v_queue_full ? local_vecs : local_pointer + 1;

                    for (size_t localOrthoCount = 0; localOrthoCount < localOrthoLimit; ++localOrthoCount) {
                        double d = std::inner_product(v, v + n, local_v + n * localOrthoCount, 0.0);
                        daxpy(n, -d, local_v + localOrthoCount * n, v);
                    }
                }

                alpha = norm(n, v);

                if (alpha > zero) {
                    scale(n, one / alpha, v);
                }
            }

            // At this point, beta = beta_{k+1}, alpha = alpha_{k+1}.


            //----------------------------------------------------------------
            // Construct rotation Qhat_{k,2k+1}.

            double alphahat = norm(alphabar, m_damp);
            double chat = alphabar / alphahat;
            double shat = m_damp / alphahat;

            // Use a plane rotation (Q_i) to turn B_i to R_i.

            double rhoold = rho;
            rho = norm(alphahat, beta);
            double c = alphahat / rho;
            double s = beta / rho;
            double thetanew = s * alpha;
            alphabar = c * alpha;

            // Use a plane rotation (Qbar_i) to turn R_i^T into R_i^bar.

            double rhobarold = rhobar;
            double zetaold = zeta;
            double thetabar = sbar * rho;
            double rhotemp = cbar * rho;
            rhobar = norm(cbar * rho, thetanew);
            cbar = cbar * rho / rhobar;
            sbar = thetanew / rhobar;
            zeta = cbar * zetabar;
            zetabar = -sbar * zetabar;

            // Update h, h_hat, x.

            for (size_t i = 0; i < n; i++) {
                hbar[i] = h[i] - (thetabar * rho / (rhoold * rhobarold)) * hbar[i];
                x[i] = x[i] + (zeta / (rho * rhobar)) * hbar[i];
                h[i] = v[i] - (thetanew / rho) * h[i];
            }

            // Estimate ||r||.

            // Apply rotation Qhat_{k,2k+1}.
            double betaacute = chat * betadd;
            double betacheck = -shat * betadd;

            // Apply rotation Q_{k,k+1}.
            double betahat = c * betaacute;
            betadd = -s * betaacute;

            // Apply rotation Qtilde_{k-1}.
            // betad = betad_{k-1} here.

            double thetatildeold = thetatilde;
            double rhotildeold = norm(rhodold, thetabar);
            double ctildeold = rhodold / rhotildeold;
            double stildeold = thetabar / rhotildeold;
            thetatilde = stildeold * rhobar;
            rhodold = ctildeold * rhobar;
            betad = -stildeold * betad + ctildeold*betahat;

            // betad   = betad_k here.
            // rhodold = rhod_k  here.

            tautildeold = (zetaold - thetatildeold * tautildeold) / rhotildeold;
            double taud = (zeta - thetatilde * tautildeold) / rhodold;
            d = d + betacheck * betacheck;
            m_norm_r = std::sqrt(d + sqr(betad - taud) + sqr(betadd));

            // Estimate ||A||.
            norm_a2 = norm_a2 + sqr(beta);
            m_norm_a = std::sqrt(norm_a2);
            norm_a2 = norm_a2 + sqr(alpha);

            // Estimate cond(A).
            maxrbar = std::max(maxrbar, rhobarold);
            if (m_itn > 1) {
                minrbar = std::min(minrbar, rhobarold);
            }
            m_cond_a = std::max(maxrbar, rhotemp) / std::min(minrbar, rhotemp);

            //----------------------------------------------------------------
            //Test for convergence.
            //---------------------------------------------------------------

            // Compute norms for convergence testing.
            m_norm_ar = std::abs(zetabar);
            m_norm_x = norm(n, x);

            // Now use these norms to estimate certain other quantities,
            // some of which will be small near a solution.

            test1 = m_norm_r / m_norm_b;
            test2 = m_norm_ar / (m_norm_a * m_norm_r);
            double test3 = one / m_cond_a;
            double t1 = test1 / (one + m_norm_a * m_norm_x / m_norm_b);
            double rtol = m_btol + m_atol * m_norm_a * m_norm_x / m_norm_b;

            // The following tests guard against extremely small values of
            // m_atol, m_btol or ctol.  (The user may have set any or all of
            // the parameters m_atol, m_btol, m_conlim  to 0.)
            // The effect is equivalent to the normAl tests using
            // m_atol = m_eps,  m_btol = m_eps,  m_conlim = 1/m_eps.

            if (m_itn >= m_itnlim) m_istop = 7;
            if (one + test3 <= one) m_istop = 6;
            if (one + test2 <= one) m_istop = 5;
            if (one + t1 <= one) m_istop = 4;

            // Allow for tolerances set by the user.

            if ( test3   <= ctol ) m_istop = 3;
            if ( test2   <= m_atol ) m_istop = 2;
            if ( test1   <= rtol ) m_istop = 1;

            //----------------------------------------------------------------
            // See if it is time to print something.
            //----------------------------------------------------------------
            if ( m_nout ) {
            bool prnt = false;
            if ( n<=40 ) prnt = true;
            if ( m_itn <= 10 ) prnt = true;
            if ( m_itn >= m_itnlim-10 ) prnt = true;
            if ( (m_itn % 10)  ==  0 ) prnt = true;
            if ( test3 <=  1.1*ctol ) prnt = true;
            if ( test2 <=  1.1*m_atol ) prnt = true;
            if ( test1 <=  1.1*rtol ) prnt = true;
            if ( m_istop!=0 ) prnt = true;

            if ( prnt ) { // Print a line for this iteration
            if ( pcount >= pfreq ) { // Print a heading first
            pcount = 0;
            if ( m_damped )
                {
                (*m_nout) << "   m_Itn       x(1)           norm rbar    Abar'rbar"
                " Compatible    LS    norm Abar cond Abar\n";
                } else {
                (*m_nout) << "   m_Itn       x(1)            norm r         A'r   "
                " Compatible    LS      norm A    cond A\n";
            }
            }
            pcount = pcount + 1;
            (*m_nout)
            << m_itn << ", " << x[0] << ", " <<m_norm_r << ", " << m_norm_ar << ", " << test1 << ", " << test2
            << ", " << m_norm_a << ", " << m_cond_a << std::endl;
            }
            }

        } while ( m_istop == 0);

        termination_print_out();
    }

private:
    void termination_print_out()
    {
        if (m_damped && m_istop == 2) {
            m_istop=3;
        }

        if (m_nout) {
            (*m_nout)
                << " Exit  LSMR. istop   = " << m_istop  << ", itn     = " << m_itn     << std::endl
                << " Exit  LSMR. norm_a  = " << m_norm_a << ", cond_a  = " << m_cond_a  << std::endl
                << " Exit  LSMR. norm_b  = " << m_norm_b << ", norm_x  = " << m_norm_x  << std::endl
                << " Exit  LSMR. norm_r  = " << m_norm_r << ", norm_ar = " << m_norm_ar << std::endl
                << " Exit  LSMR. " << stopping_reason_message() << std::endl;
        }
    }

    double m_norm_a;
    double m_cond_a;
    double m_norm_b;
    double m_norm_r;
    double m_norm_ar;
    double m_norm_x;
    double m_dxmax;

    double m_atol;
    double m_btol;
    double m_conlim;

    double m_eps;
    double m_damp;
    bool m_damped;

    size_t m_itnlim;
    size_t m_itn;

    size_t m_istop;

    size_t m_maxdx;
    size_t m_local_size;
    std::ostream* m_nout;
};

template <typename TMatrix>
class LsmrSymmetricSolver : public LsmrSolverBase
{
private:
    Ref<const TMatrix> m_a;

public:
    LsmrSymmetricSolver(Ref<const TMatrix> a) : m_a(a) { }

    ~LsmrSymmetricSolver() { }

    void a_prod_1(Ref<const Vector> x, Ref<Vector> y) const override
    {
        y += m_a.selfadjointView<Eigen::Upper>() * x;
    }

    void a_prod_2(Ref<Vector> x, Ref<const Vector> y) const override
    {
        x += m_a.selfadjointView<Eigen::Upper>() * y;
    }

    void analyze_pattern(Ref<const TMatrix> value) override
    { }

    void set_matrix(Ref<const TMatrix> value) override
    {
        m_a = value;
    }
};

using LsmrDenseSolver = LsmrSymmetricSolver<Matrix>;
using SolverLSMR = LsmrSymmetricSolver<Sparse>;

} // namespace EQlib