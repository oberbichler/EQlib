#pragma once

#include "Define.h"
#include "Timer.h"

#include "tbb/blocked_range.h"
#include "tbb/parallel_reduce.h"

namespace EQlib {

struct Assemble
{
    double m_computation_time;

    py::dict& m_options;

    Vector m_lhs_values;
    Vector m_rhs_values;

    Map<Sparse> m_lhs;
    Map<Vector> m_rhs;

    Assemble(Sparse& lhs, Vector& rhs, py::dict& options)
    : m_lhs(lhs.rows(), lhs.cols(), lhs.nonZeros(), lhs.outerIndexPtr(),
        lhs.innerIndexPtr(), lhs.valuePtr())
    , m_rhs(rhs.data(), rhs.size())
    , m_options(options)
    , m_computation_time(0)
    {
        // set lhs and rhs to zero
        Map<Vector>(m_lhs.valuePtr(), m_lhs.nonZeros()).setZero();
        m_rhs.setZero();
    }

    Assemble(Assemble& s, tbb::split)
    : m_lhs_values(Vector::Zero(s.m_lhs.nonZeros()))
    , m_rhs_values(Vector::Zero(s.m_rhs.size()))
    , m_lhs(s.m_lhs.rows(), s.m_lhs.cols(), s.m_lhs.nonZeros(),
        s.m_lhs.outerIndexPtr(), s.m_lhs.innerIndexPtr(),
        m_lhs_values.data())
    , m_rhs(m_rhs_values.data(), s.m_rhs.size())
    , m_options(s.m_options)
    , m_computation_time(0)
    { }

    template <typename TRange>
    void operator()(const TRange& range)
    {
        Timer timer;

        // compute and add local lhs and rhs

        for (auto it = range.begin(); it != range.end(); ++it) {
            const auto& [element, dof_indices] = *it;

            const auto [local_lhs, local_rhs] = element->compute(m_options);

            const size_t nb_dofs = dof_indices.size();

            for (size_t row = 0; row < nb_dofs; row++) {
                const auto row_index = dof_indices[row];

                if (row_index.global >= m_rhs.size()) {
                    continue;
                }

                m_rhs(row_index.global) += local_rhs(row_index.local);

                for (size_t col = row; col < nb_dofs; col++) {
                    const auto col_index = dof_indices[col];

                    if (col_index.global >= m_rhs.size()) {
                        continue;
                    }

                    m_lhs.coeffRef(row_index.global, col_index.global) +=
                        local_lhs(row_index.local, col_index.local);
                }
            }
        }

        m_computation_time = timer.ellapsed();
    }

    void join(Assemble& rhs)
    {
        Map<Vector>(m_lhs.valuePtr(), m_lhs.nonZeros()) += rhs.m_lhs_values;
        Map<Vector>(m_rhs.data(), m_rhs.size()) += rhs.m_rhs_values;
        m_computation_time += rhs.m_computation_time;
    }

    template <typename TContainer>
    static void serial(TContainer& element_indices, Sparse& lhs, Vector& rhs,
        py::dict& options)
    {
        Assemble result(lhs, rhs, options);

        result(element_indices);
    }

    template <typename TContainer>
    static void parallel(TContainer& element_indices, Sparse& lhs,
        Vector& rhs, py::dict& options)
    {
        py::gil_scoped_release release;

        Assemble result(lhs, rhs, options);

        tbb::parallel_reduce(tbb::blocked_range<TContainer::iterator>(
            element_indices.begin(), element_indices.end()), result);
    }
};

} // namespace EQlib