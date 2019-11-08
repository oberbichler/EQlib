#pragma once

#include "Define.h"

namespace eqlib {

class ProblemData
{
private:    // variables
    index m_n;
    index m_m;
    index m_nb_nonzeros_dg;
    index m_nb_nonzeros_hl;
    std::vector<double> m_values;

public:     // variables
    double m_computation_time;
    double m_assemble_time;
    std::vector<double> m_buffer;

public:     // constructor
    ProblemData() : m_computation_time(0), m_assemble_time(0)
    {
    }

public:     // methods
    double& computation_time()
    {
        return m_computation_time;
    }

    double& assemble_time()
    {
        return m_assemble_time;
    }

    void set_zero()
    {
        std::fill(m_values.begin(), m_values.end(), 0);
        std::fill(m_buffer.begin(), m_buffer.end(), 0);
        m_computation_time = 0.0;
        m_assemble_time = 0.0;
    }

    void resize(const index n, const index m, const index nb_nonzeros_dg, const index nb_nonzeros_hl, const index max_element_n, const index max_element_m)
    {
        m_n = n;
        m_m = m;
        m_nb_nonzeros_dg = nb_nonzeros_dg;
        m_nb_nonzeros_hl = nb_nonzeros_hl;

        const index nb_entries = 1 + m + n + nb_nonzeros_dg + nb_nonzeros_hl;

        m_values.resize(nb_entries);

        m_buffer.resize(
            std::max(max_element_m, max_element_n) +
            max_element_m * max_element_n +
            max_element_n * max_element_n);

        set_zero();
    }

    double& f() noexcept
    {
        return m_values[0];
    }

    double f() const noexcept
    {
        return m_values[0];
    }

    double& g(const index i) noexcept
    {
        assert(0 <= i && i < m_m);
        return m_values[1 + i];
    }

    double g(const index i) const noexcept
    {
        assert(0 <= i && i < m_m);
        return m_values[1 + i];
    }

    Map<Vector> g() noexcept
    {
        return Map<Vector>(m_values.data() + 1, m_m);
    }

    Map<const Vector> g() const noexcept
    {
        return Map<const Vector>(m_values.data() + 1, m_m);
    }

    double& df(const index i) noexcept
    {
        assert(0 <= i && i < m_n);
        return m_values[1 + m_m + i];
    }

    double df(const index i) const noexcept
    {
        assert(0 <= i && i < m_n);
        return m_values[1 + m_m + i];
    }

    Map<Vector> df() noexcept
    {
        return Map<Vector>(m_values.data() + 1 + m_m, m_n);
    }

    Map<const Vector> df() const noexcept
    {
        return Map<const Vector>(m_values.data() + 1 + m_m, m_n);
    }

    double& dg(const index i) noexcept
    {
        assert(0 <= i && i < m_nb_nonzeros_dg);
        return m_values[1 + m_m + m_n + i];
    }

    double dg(const index i) const noexcept
    {
        assert(0 <= i && i < m_nb_nonzeros_dg);
        return m_values[1 + m_m + m_n + i];
    }

    Map<Vector> dg() noexcept
    {
        return Map<Vector>(m_values.data() + 1 + m_m + m_n, m_nb_nonzeros_dg);
    }

    Map<const Vector> dg() const noexcept
    {
        return Map<const Vector>(m_values.data() + 1 + m_m + m_n, m_nb_nonzeros_dg);
    }

    const double* dg_ptr() const noexcept
    {
        return m_values.data() + 1 + m_m + m_n;
    }

    double& hm(const index i) noexcept
    {
        assert(0 <= i && i < m_nb_nonzeros_hl);
        return m_values[1 + m_m + m_n + m_nb_nonzeros_dg + i];
    }

    double hm(const index i) const noexcept
    {
        assert(0 <= i && i < m_nb_nonzeros_hl);
        return m_values[1 + m_m + m_n + m_nb_nonzeros_dg + i];
    }

    const double* const hl_ptr() const noexcept
    {
        return m_values.data() + 1 + m_m + m_n + m_nb_nonzeros_dg;
    }

    Map<Vector> hm() noexcept
    {
        return Map<Vector>(m_values.data() + 1 + m_m + m_n + m_nb_nonzeros_dg, m_nb_nonzeros_hl);
    }

    Map<const Vector> hm() const noexcept
    {
        return Map<const Vector>(m_values.data() + 1 + m_m + m_n + m_nb_nonzeros_dg, m_nb_nonzeros_hl);
    }

    const std::vector<double> values() const
    {
        return m_values;
    }

    double* values_ptr()
    {
        return m_values.data();
    }

    const double* values_ptr() const
    {
        return m_values.data();
    }

    ProblemData& operator+=(const ProblemData& rhs)
    {
        assert(length(m_values) == length(rhs.m_values));

        for (index i = 0; i < length(m_values); i++) {
            m_values[i] += rhs.m_values[i];
        }

        m_computation_time += rhs.m_computation_time;
        m_assemble_time += rhs.m_assemble_time;

        return *this;
    }
};

} // namespace eqlib