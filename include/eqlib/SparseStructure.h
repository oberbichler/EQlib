#pragma once

#include "Define.h"

#include <sparsehash/dense_hash_map>

#include <vector>

namespace eqlib {

template <typename TScalar = double, typename TIndex = int, bool TRowMajor = false, bool TIndexMap = true>
class SparseStructure
{
private:    // variables
    TIndex m_rows;
    TIndex m_cols;
    std::vector<TIndex> m_ia;
    std::vector<TIndex> m_ja;
    std::vector<google::dense_hash_map<TIndex, index>> m_indices;

public:     // methods
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
        const index size = rows() * cols();

        if (size == 0) {
            return 0;
        }

        return (double)nb_nonzeros() / size;
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

    index get_index(const index row, const index col) const
    {
        assert(0 <= row && row < rows());
        assert(0 <= col && col < cols());

        const auto i = static_cast<TIndex>(TRowMajor ? row : col);
        const auto j = static_cast<TIndex>(TRowMajor ? col : row);

        if (TIndexMap) {
            const auto it = m_indices[i].find(j);

            assert(it != m_indices[i].end());

            return it->second;
        } else {
            const auto lower = m_ja.begin() + m_ia[i];
            const auto upper = m_ja.begin() + m_ia[i + 1];

            const auto it = std::lower_bound(lower, upper, j);

            if (*it != j || it == upper) {
                assert(false);
                return -1;
            }

            const index value_index = std::distance(m_ja.begin(), it);

            assert(value_index < nb_nonzeros());

            return value_index;
        }
    }

    template <typename TPattern>
    void set(const index rows, const index cols, const TPattern& pattern) noexcept
    {
        assert(length(pattern) == (TRowMajor ? rows : cols));

        m_rows = static_cast<TIndex>(rows);
        m_cols = static_cast<TIndex>(cols);

        m_ia.resize((TRowMajor ? rows : cols) + 1);

        m_ia[0] = 0;

        for (TIndex i = 0; i < (TRowMajor ? rows : cols); i++) {
            const TIndex n = static_cast<TIndex>(pattern[i].size());

            m_ia[i + 1] = m_ia[i] + n;
        }

        m_ja.resize(nb_nonzeros());

        auto ja_it = m_ja.begin();

        for (TIndex i = 0; i < (TRowMajor ? rows : cols); i++) {
            const TIndex n = static_cast<TIndex>(pattern[i].size());

            for (const auto j : pattern[i]) {
                assert(j < (TRowMajor ? m_cols : m_rows));
                *ja_it++ = static_cast<TIndex>(j);
            }
        }

        if (TIndexMap) {
            m_indices.resize(TRowMajor ? rows : cols);

            for (index i = 0; i < (TRowMajor ? rows : cols); i++) {
                m_indices[i].set_empty_key(-1);
                m_indices[i].resize(m_ia[i + 1] - m_ia[i]);
                for (TIndex k = m_ia[i]; k < m_ia[i + 1]; k++) {
                    const TIndex j = m_ja[k];
                    m_indices[i][j] = k;
                }
            }
        }
    }

    /*
    * https://github.com/scipy/scipy/blob/3b36a574dc657d1ca116f6e230be694f3de31afc/scipy/sparse/sparsetools/csr.h#L380
    */
    void convert_from(SparseStructure<TScalar, TIndex, !TRowMajor> other, std::vector<TScalar>& values)
    {
        const auto nnz = other.nb_nonzeros();

        const auto n = TRowMajor ? other.cols() : other.rows();
        const auto m = TRowMajor ? other.rows() : other.cols();

        m_rows = other.rows();
        m_cols = other.cols();

        m_ia.resize(m + 1);
        m_ja.resize(nnz);

        std::fill(m_ia.begin(), m_ia.end(), 0);

        for (TIndex n = 0; n < nnz; n++) {
            m_ia[other.ja(n)]++;
        }

        TIndex cumsum = 0;

        for (TIndex j = 0; j < m; j++) {
            const auto temp  = m_ia[j];
            m_ia[j] = cumsum;
            cumsum += temp;
        }

        m_ia[m] = nnz;

        std::vector<TScalar> a_values {values};

        for (TIndex i = 0; i < n; i++){
            for(TIndex k = other.ia(i); k < other.ia(i + 1); k++){
                const auto j = other.ja(k);
                const auto dest = m_ia[j];

                m_ja[dest] = i;
                values[dest] = a_values[k];

                m_ia[j]++;
            }
        }

        TIndex last = 0;

        for (TIndex j = 0; j <= m; j++){
            const auto temp = m_ia[j];
            m_ia[j] = last;
            last = temp;
        }
    }

    void for_each(std::function<void(TIndex, TIndex)> action) const
    {
        for (TIndex col = 0; col < m_cols; col++) {
            for (TIndex i = m_ia[col]; i < m_ia[col + 1]; i++) {
                const TIndex row = m_ja[i];
                action(row, col);
            }
        }
    }
};

} // namespace eqlib