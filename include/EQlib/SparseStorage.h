#pragma once

#include "Define.h"

#include <vector>

namespace EQlib {

template <typename TScalar = double, typename TIndex = int, bool TRowMajor = false>
class SparseStorage
{
private:    // variables
    TIndex m_rows;
    TIndex m_cols;
    std::vector<TIndex> m_ia;
    std::vector<TIndex> m_ja;

public:     // methods
    TIndex rows() const noexcept
    {
        return m_rows;
    }

    TIndex cols() const noexcept
    {
        return m_cols;
    }

    TIndex nnz() const noexcept
    {
        return m_ia.back();
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


    template <typename TContainer>
    TScalar& coeff_ref(TContainer& values, const index row, const index col) noexcept
    {
        assert(length(values) == nnz());
        assert(row < rows());
        assert(col < cols());

        static TScalar dummy;

        TIndex i = static_cast<TIndex>(TRowMajor ? row : col);
        TIndex j = static_cast<TIndex>(TRowMajor ? col : row);

        const auto lower = m_ja.begin() + m_ia[i];
        const auto upper = m_ja.begin() + m_ia[i + 1];

        const auto it = std::lower_bound(lower, upper, j);

        if (*it != j || it == upper) {
            return dummy;
        }

        const auto index = std::distance(m_ja.begin(), it);

        return values[index];
    }

    template <typename TContainer>
    TScalar coeff(const TContainer& values, const index row, const index col) const noexcept
    {
        assert(length(values) == nnz());
        assert(row < rows());
        assert(col < cols());

        TIndex i = static_cast<TIndex>(TRowMajor ? row : col);
        TIndex j = static_cast<TIndex>(TRowMajor ? col : row);

        const auto lower = m_ja.begin() + m_ia[i];
        const auto upper = m_ja.begin() + m_ia[i + 1];

        const auto it = std::lower_bound(lower, upper, j);

        if (*it != j || it == upper) {
            return 0;
        }

        const auto index = std::distance(m_ja.begin(), it);

        return values[index];
    }

    template <typename TPattern>
    void set(const index rows, const index cols, const TPattern& pattern) noexcept
    {
        assert(length(pattern) == (TRowMajor ? rows : cols));

        m_rows = static_cast<TIndex>(rows);
        m_cols = static_cast<TIndex>(cols);

        m_ia.resize((TRowMajor ? rows : cols) + 1);

        m_ia[0] = 0;

        for (TIndex i = 0; i < cols; i++) {
            const TIndex n = static_cast<TIndex>(pattern[i].size());

            m_ia[i + 1] = m_ia[i] + n;
        }

        m_ja.resize(nnz());

        auto ja_it = m_ja.begin();

        for (TIndex i = 0; i < cols; i++) {
            const TIndex n = static_cast<TIndex>(pattern[i].size());

            for (const auto j : pattern[i]) {
                assert(j < (TRowMajor ? m_cols : m_rows));
                *ja_it++ = static_cast<TIndex>(j);
            }
            
            std::sort(ja_it - n, ja_it);
        }
    }

    /*
    * https://github.com/scipy/scipy/blob/3b36a574dc657d1ca116f6e230be694f3de31afc/scipy/sparse/sparsetools/csr.h#L380
    */
    void convert_from(SparseStorage<TScalar, TIndex, !TRowMajor> other, std::vector<TScalar>& values)
    {
        const auto nnz = other.nnz();

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
};

} // namespace EQlib