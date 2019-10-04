#pragma once

#include "Define.h"

namespace EQlib {

class ProblemData
{
private:    // variables
    Vector m_data;
    index m_n;
    index m_m;
    index m_nnz_dg;
    index m_nnz_hl;

public:     // methods
    void set_zero()
    {
        m_data.setZero();
    }

    void set_zero(const index n, const index m, const index nnz_dg, const index nnz_hl)
    {
        m_n = n;
        m_m = m;
        m_nnz_dg = nnz_dg;
        m_nnz_hl = nnz_hl;

        const index nb_entries = 1 + m + n + nnz_dg + nnz_hl;

        m_data.resize(nb_entries);

        m_data.setZero();
    }

    double& f() noexcept
    {
        return m_data[0];
    }

    double f() const noexcept
    {
        return m_data[0];
    }

    double& g(const index i) noexcept
    {
        return m_data[1 + i];
    }

    double g(const index i) const noexcept
    {
        return m_data[1 + i];
    }

    Map<Vector> g() noexcept
    {
        return Map<Vector>(m_data.data() + 1, m_m);
    }

    Map<const Vector> g() const noexcept
    {
        return Map<const Vector>(m_data.data() + 1, m_m);
    }

    double& df(const index i) noexcept
    {
        return m_data[1 + m_m + i];
    }

    double df(const index i) const noexcept
    {
        return m_data[1 + m_m + i];
    }

    Map<Vector> df() noexcept
    {
        return Map<Vector>(m_data.data() + 1 + m_m, m_n);
    }

    Map<const Vector> df() const noexcept
    {
        return Map<const Vector>(m_data.data() + 1 + m_m, m_n);
    }

    double& dg(const index i) noexcept
    {
        return m_data[1 + m_m + m_n + i];
    }

    double dg(const index i) const noexcept
    {
        return m_data[1 + m_m + m_n + i];
    }

    Map<Vector> dg() noexcept
    {
        return Map<Vector>(m_data.data() + 1 + m_m + m_n, m_nnz_dg);
    }

    Map<const Vector> dg() const noexcept
    {
        return Map<const Vector>(m_data.data() + 1 + m_m + m_n, m_nnz_dg);
    }

    const double* dg_ptr() const noexcept
    {
        return m_data.data() + 1 + m_m + m_n;
    }

    double& hl(const index i) noexcept
    {
        return m_data[1 + m_m + m_n + m_nnz_dg + i];
    }

    double hl(const index i) const noexcept
    {
        return m_data[1 + m_m + m_n + m_nnz_dg + i];
    }

    const double* const hl_ptr() const noexcept
    {
        return m_data.data() + 1 + m_m + m_n + m_nnz_dg;
    }

    Map<Vector> hl() noexcept
    {
        return Map<Vector>(m_data.data() + 1 + m_m + m_n + m_nnz_dg, m_nnz_hl);
    }

    Map<const Vector> hl() const noexcept
    {
        return Map<const Vector>(m_data.data() + 1 + m_m + m_n + m_nnz_dg, m_nnz_hl);
    }

    ProblemData& operator+=(const ProblemData& rhs)
    {
        m_data += rhs.m_data;
        return *this;
    }
};

} // namespace EQlib