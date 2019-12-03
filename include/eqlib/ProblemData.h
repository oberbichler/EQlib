#pragma once

#include "Define.h"

namespace eqlib {

class ProblemData {
private: // variables
    index m_n;
    index m_m;
    index m_nb_nonzeros_dg;
    index m_nb_nonzeros_hm;
    Vector m_values;

public: // variables
    double m_computation_time;
    double m_assemble_time;
    Vector m_buffer;

public: // constructor
    ProblemData()
        : m_computation_time(0)
        , m_assemble_time(0)
    {
    }

public: // methods
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
        m_values.setZero();
        m_buffer.setZero();
        m_computation_time = 0.0;
        m_assemble_time = 0.0;
    }

    void resize(
        const index n,
        const index m,
        const index nb_nonzeros_dg,
        const index nb_nonzeros_hm,
        const index max_element_n,
        const index max_element_m)
    {
        m_n = n;
        m_m = m;
        m_nb_nonzeros_dg = nb_nonzeros_dg;
        m_nb_nonzeros_hm = nb_nonzeros_hm;

        const index nb_entries = 1 + m + n + nb_nonzeros_dg + nb_nonzeros_hm;

        m_values.resize(nb_entries);

        m_buffer.resize(
            std::max(index{1}, max_element_m) * max_element_n + std::max(index{1}, max_element_m) * max_element_n * max_element_n);

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

    Ref<Vector> g() noexcept
    {
        return m_values.segment(1, m_m);
    }

    Ref<const Vector> g() const noexcept
    {
        return m_values.segment(1, m_m);
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

    Ref<Vector> df() noexcept
    {
        return m_values.segment(1 + m_m, m_n);
    }

    Ref<const Vector> df() const noexcept
    {
        return m_values.segment(1 + m_m, m_n);
    }

    double& dg_value(const index i) noexcept
    {
        assert(0 <= i && i < m_nb_nonzeros_dg);
        return m_values[1 + m_m + m_n + i];
    }

    double dg_value(const index i) const noexcept
    {
        assert(0 <= i && i < m_nb_nonzeros_dg);
        return m_values[1 + m_m + m_n + i];
    }

    Ref<Vector> dg_values() noexcept
    {
        return m_values.segment(1 + m_m + m_n, m_nb_nonzeros_dg);
    }

    Ref<const Vector> dg_values() const noexcept
    {
        return m_values.segment(1 + m_m + m_n, m_nb_nonzeros_dg);
    }

    const double* dg_ptr() const noexcept
    {
        return m_values.data() + 1 + m_m + m_n;
    }

    double& hm_value(const index i) noexcept
    {
        assert(0 <= i && i < m_nb_nonzeros_hm);
        return m_values[1 + m_m + m_n + m_nb_nonzeros_dg + i];
    }

    double hm_value(const index i) const noexcept
    {
        assert(0 <= i && i < m_nb_nonzeros_hm);
        return m_values[1 + m_m + m_n + m_nb_nonzeros_dg + i];
    }

    const double* const hm_ptr() const noexcept
    {
        return m_values.data() + 1 + m_m + m_n + m_nb_nonzeros_dg;
    }

    Ref<Vector> hm_values() noexcept
    {
        return m_values.segment(1 + m_m + m_n + m_nb_nonzeros_dg, m_nb_nonzeros_hm);
    }

    Ref<const Vector> hm_values() const noexcept
    {
        return m_values.segment(1 + m_m + m_n + m_nb_nonzeros_dg, m_nb_nonzeros_hm);
    }

    Ref<Vector> values()
    {
        return m_values;
    }

    Ref<const Vector> values() const
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

        m_values += rhs.m_values;

        m_computation_time += rhs.m_computation_time;
        m_assemble_time += rhs.m_assemble_time;

        return *this;
    }
};

} // namespace eqlib