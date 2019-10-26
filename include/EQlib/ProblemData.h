#pragma once

#include "Define.h"

namespace EQlib {

class ProblemData
{
private:    // variables
    index m_n;
    index m_m;
    index m_nb_nonzeros_dg;
    index m_nb_nonzeros_hl;

public:     // variables
    double m_timer_allocate;
    double m_timer_compute;
    double m_timer_assemble;
    std::vector<double> m_buffer;
    std::vector<double> m_data;

public:     // constructor
    ProblemData() : m_timer_allocate(0), m_timer_compute(0), m_timer_assemble(0)
    {
    }

public:     // methods
    void set_zero()
    {
        std::fill(m_data.begin(), m_data.end(), 0);
        std::fill(m_buffer.begin(), m_buffer.end(), 0);
    }

    void resize(const index n, const index m, const index nb_nonzeros_dg, const index nb_nonzeros_hl, const index max_element_n, const index max_element_m)
    {
        m_n = n;
        m_m = m;
        m_nb_nonzeros_dg = nb_nonzeros_dg;
        m_nb_nonzeros_hl = nb_nonzeros_hl;

        const index nb_entries = 1 + m + n + nb_nonzeros_dg + nb_nonzeros_hl;

        m_data.resize(nb_entries);

        m_buffer.resize(
            std::max(max_element_m , max_element_n) +
            max_element_m * max_element_n +
            max_element_n * max_element_n);

        set_zero();
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
        return Map<Vector>(m_data.data() + 1 + m_m + m_n, m_nb_nonzeros_dg);
    }

    Map<const Vector> dg() const noexcept
    {
        return Map<const Vector>(m_data.data() + 1 + m_m + m_n, m_nb_nonzeros_dg);
    }

    const double* dg_ptr() const noexcept
    {
        return m_data.data() + 1 + m_m + m_n;
    }

    double& hl(const index i) noexcept
    {
        return m_data[1 + m_m + m_n + m_nb_nonzeros_dg + i];
    }

    double hl(const index i) const noexcept
    {
        return m_data[1 + m_m + m_n + m_nb_nonzeros_dg + i];
    }

    const double* const hl_ptr() const noexcept
    {
        return m_data.data() + 1 + m_m + m_n + m_nb_nonzeros_dg;
    }

    Map<Vector> hl() noexcept
    {
        return Map<Vector>(m_data.data() + 1 + m_m + m_n + m_nb_nonzeros_dg, m_nb_nonzeros_hl);
    }

    Map<const Vector> hl() const noexcept
    {
        return Map<const Vector>(m_data.data() + 1 + m_m + m_n + m_nb_nonzeros_dg, m_nb_nonzeros_hl);
    }

    ProblemData& operator+=(const ProblemData& rhs)
    {
        assert(length(m_data) == length(rhs.m_data));

        for (index i = 0; i < length(m_data); i++) {
            m_data[i] += rhs.m_data[i];
        }

        m_timer_allocate += rhs.m_timer_allocate;
        m_timer_compute += rhs.m_timer_compute;
        m_timer_assemble += rhs.m_timer_assemble;

        return *this;
    }
};

} // namespace EQlib