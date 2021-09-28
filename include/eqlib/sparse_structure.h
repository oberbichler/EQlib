#pragma once

#include "common.h"

#include <set>
#include <vector>

namespace eqlib {

template <typename TScalar = double, typename TIndex = int, bool TRowMajor = true>
struct SparseStructure {
    TIndex m_rows;
    TIndex m_cols;
    std::vector<TIndex> m_ia;
    std::vector<TIndex> m_ja;

    SparseStructure()
    {
    }

    SparseStructure(const TIndex rows, const TIndex cols, const std::vector<TIndex>& ia, const std::vector<TIndex>& ja)
        : m_rows(rows)
        , m_cols(cols)
        , m_ia(ia)
        , m_ja(ja)
    {
        const TIndex size_i = TRowMajor ? rows : cols;
        const TIndex size_j = TRowMajor ? cols : rows;

        if (len(ia) != size_i + 1) {
            throw std::invalid_argument("Vector ia has an invalid size");
        }

        if (ja.size() > 0) {
            const TIndex max_j = Map<const Eigen::Matrix<TIndex, 1, Eigen::Dynamic>>(ja.data(), ja.size()).maxCoeff();

            if (max_j >= size_j) {
                throw std::invalid_argument("Vector ja has invalid entries");
            }
        }
    }

    TIndex rows() const noexcept
    {
        return m_rows;
    }

    TIndex cols() const noexcept
    {
        return m_cols;
    }

    TIndex nb_nonzeros() const noexcept
    {
        return m_ia.back();
    }

    double density() const noexcept
    {
        if (rows() == 0 || cols() == 0) {
            return 0;
        }

        return (double)nb_nonzeros() / rows() / cols();
    }

    const std::vector<TIndex>& ia() const noexcept
    {
        return m_ia;
    }

    const TIndex& ia(TIndex index) const
    {
        return m_ia.at(index);
    }

    const std::vector<TIndex>& ja() const noexcept
    {
        return m_ja;
    }

    const TIndex& ja(TIndex index) const
    {
        return m_ja.at(index);
    }

    index get_first_index(const index i) const
    {
        return static_cast<index>(m_ia[i]);
    }

    index get_index_bounded(const index j, const index lo, const index hi) const
    {
        const auto ja_begin = m_ja.begin();

        const auto lower = std::next(ja_begin, lo);
        const auto upper = std::next(ja_begin, hi);

        const auto it = std::lower_bound(lower, upper, j);

        if (*it != j || it == upper) {
            return -1;
        }

        const index value_index = std::distance(ja_begin, it);

        assert(value_index < nb_nonzeros());

        return value_index;
    }

    index get_index(const index row, const index col) const
    {
        assert(0 <= row && row < rows());
        assert(0 <= col && col < cols());

        const auto i = static_cast<TIndex>(TRowMajor ? row : col);
        const auto j = static_cast<TIndex>(TRowMajor ? col : row);

        const auto ja_begin = m_ja.begin();

        const auto lower = std::next(ja_begin, m_ia[i]);
        const auto upper = std::next(ja_begin, m_ia[i + 1]);

        const auto it = std::lower_bound(lower, upper, j);

        if (*it != j || it == upper) {
            return -1;
        }

        const index value_index = std::distance(ja_begin, it);

        assert(value_index < nb_nonzeros());

        return value_index;
    }

    template <typename TPattern>
    static auto from_pattern(const index rows, const index cols, const TPattern& pattern)
    {
        const index size_i = TRowMajor ? rows : cols;
        const index size_j = TRowMajor ? cols : rows;

        assert(len(pattern) == size_i);

        std::vector<TIndex> ia(size_i + 1);

        ia[0] = 0;

        for (TIndex i = 0; i < size_i; i++) {
            const TIndex n = static_cast<TIndex>(pattern[i].size());

            ia[i + 1] = ia[i] + n;
        }

        const auto nb_nonzeros = ia.back();

        std::vector<TIndex> ja(nb_nonzeros);

        auto ja_it = ja.begin();

        for (TIndex i = 0; i < size_i; i++) {
            const TIndex n = static_cast<TIndex>(pattern[i].size());

            for (const auto j : pattern[i]) {
                assert(j < size_j);
                *ja_it++ = static_cast<TIndex>(j);
            }

            std::sort(ja_it - n, ja_it);
        }

        return SparseStructure<TScalar, TIndex, TRowMajor>(static_cast<TIndex>(rows), static_cast<TIndex>(cols), ia, ja);
    }

    std::pair<SparseStructure, std::vector<index>> to_general() const
    {
        assert(m_rows == m_cols);

        const index n = m_rows;

        std::vector<std::vector<TIndex>> pattern(n);
        std::vector<std::vector<TIndex>> indices(n);

        for (TIndex i = 0; i < n; i++) {
            for (TIndex k = m_ia[i]; k < m_ia[i + 1]; k++) {
                const TIndex j = m_ja[k];

                pattern[i].push_back(j);
                indices[i].push_back(k);

                if (i == j) {
                    continue;
                }

                pattern[j].push_back(i);
                indices[j].push_back(k);
            }
        }

        SparseStructure result = from_pattern(m_rows, m_cols, pattern);

        std::vector<index> value_indices;
        value_indices.reserve(nb_nonzeros() * 2 - n);

        for (const auto& indices_row : indices) {
            value_indices.insert(value_indices.end(), indices_row.begin(), indices_row.end());
        }

        return { result, value_indices };
    }

    std::pair<SparseStructure, Vector> to_general(Ref<const Vector> values) const
    {
        assert(m_rows == m_cols);

        const index n = m_rows;

        std::vector<std::vector<TIndex>> pattern(n);
        std::vector<std::vector<TIndex>> indices(n);

        for (TIndex i = 0; i < n; i++) {
            for (TIndex k = m_ia[i]; k < m_ia[i + 1]; k++) {
                const TIndex j = m_ja[k];

                pattern[i].push_back(j);
                indices[i].push_back(k);

                if (i == j) {
                    continue;
                }

                pattern[j].push_back(i);
                indices[j].push_back(k);
            }
        }

        SparseStructure result = from_pattern(m_rows, m_cols, pattern);

        Vector new_values(nb_nonzeros() * 2 - n);
        index i = 0;

        for (const auto& indices_row : indices) {
            for (const auto& index : indices_row) {
                new_values(i++) = values(index);
            }
        }

        return { result, new_values };
    }
};

using CsrStructure = SparseStructure<double, int, true>;

} // namespace eqlib