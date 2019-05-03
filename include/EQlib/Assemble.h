#pragma once

#include "Define.h"
#include "Timer.h"

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/iterators.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_reduce.h>

namespace EQlib {

struct Assemble
{
    static inline tbb::enumerable_thread_specific<Vector> m_h_values;
    static inline tbb::enumerable_thread_specific<Vector> m_g_values;

    double m_f;
    Map<Vector> m_g;
    Map<Sparse> m_h;

    static inline Map<Vector> create_g(Assemble& s) {
        auto& g_values = m_g_values.local();

        g_values.resize(s.m_g.size());
        g_values.setZero();

        return Map<Vector>(g_values.data(), s.m_g.size());
    }

    static inline Map<Sparse> create_h(Assemble& s) {
        auto& h_values = m_h_values.local();

        h_values.resize(s.m_h.nonZeros());
        h_values.setZero();

        return Map<Sparse>(s.m_h.rows(), s.m_h.cols(), s.m_h.nonZeros(),
            s.m_h.outerIndexPtr(), s.m_h.innerIndexPtr(),
            h_values.data());
    }

    Assemble(Vector& g, Sparse& h)
    : m_f(0.0)
    , m_g(g.data(), g.size())
    , m_h(h.rows(), h.cols(), h.nonZeros(), h.outerIndexPtr(),
        h.innerIndexPtr(), h.valuePtr())
    {
        // set h and g to zero
        Map<Vector>(m_h.valuePtr(), m_h.nonZeros()).setZero();
        m_g.setZero();
    }

    Assemble(Assemble& s, tbb::split)
    : m_f(0.0)
    , m_g(create_g(s))
    , m_h(create_h(s))
    { }

    template <typename TRange>
    void operator()(const TRange& range)
    {
        // compute and add local h and g

        for (auto it = range.begin(); it != range.end(); ++it) {
            const auto& [element, dof_indices] = *it;

            const auto [local_f, local_g, local_h] = element->compute();

            const size_t nb_dofs = dof_indices.size();

            m_f += local_f;

            for (size_t row = 0; row < nb_dofs; row++) {
                const auto row_index = dof_indices[row];

                if (row_index.global >= m_g.size()) {
                    continue;
                }

                m_g(row_index.global) += local_g(row_index.local);

                for (size_t col = row; col < nb_dofs; col++) {
                    const auto col_index = dof_indices[col];

                    if (col_index.global >= m_g.size()) {
                        continue;
                    }

                    m_h.coeffRef(row_index.global, col_index.global) +=
                        local_h(row_index.local, col_index.local);
                }
            }
        }
    }

    void join(Assemble& rhs)
    {
        m_f += rhs.m_f;
        Map<Vector>(m_g.data(), m_g.size()) += 
            Map<Vector>(rhs.m_g.data(), rhs.m_g.size());
        Map<Vector>(m_h.valuePtr(), m_h.nonZeros()) +=
            Map<Vector>(rhs.m_h.valuePtr(), rhs.m_h.nonZeros());
    }

    template <typename TElements, typename TIndices>
    static double run(
        int nb_threads,
        TElements& elements,
        TIndices& indices,
        Vector& g,
        Sparse& h)
    {
        py::gil_scoped_release release;

        tbb::task_scheduler_init init(nb_threads > 0 ? nb_threads :
            tbb::task_scheduler_init::automatic);

        using Iterator = tbb::zip_iterator<TElements::iterator,
            TIndices::iterator>;

        auto begin = tbb::make_zip_iterator(elements.begin(), indices.begin());
        auto end = tbb::make_zip_iterator(elements.end(), indices.end());

        Assemble result(g, h);

        tbb::parallel_reduce(tbb::blocked_range<Iterator>(begin, end), result);

        return result.m_f;
    }
};

} // namespace EQlib