#pragma once

#include "Define.h"

#include <vector>

namespace eqlib {

template <typename TScalar = double, typename TIndex = int, bool TRowMajor = false, bool TIndexMap = true>
class SparseStructure
{
private:    // types
    using Type = SparseStructure<TScalar, TIndex, TRowMajor, TIndexMap>;

private:    // variables
    TIndex m_rows;
    TIndex m_cols;
    std::vector<TIndex> m_ia;
    std::vector<TIndex> m_ja;
    std::vector<DenseMap<TIndex, index>> m_indices;

public: // constructors
    SparseStructure()
    {
    }

    SparseStructure(const TIndex rows, const TIndex cols, const std::vector<TIndex>& ia, const std::vector<TIndex>& ja)
        : m_rows(rows)
        , m_cols(cols)
        , m_ia(ia)
        , m_ja(ja)
    {
        const index size_i = TRowMajor ? rows : cols;
        const index size_j = TRowMajor ? cols : rows;

        if (TIndexMap) {
            m_indices.resize(size_i);

            for (index i = 0; i < size_i; i++) {
                m_indices[i].set_empty_key(-1);
                m_indices[i].resize(m_ia[i + 1] - m_ia[i]);
                for (TIndex k = m_ia[i]; k < m_ia[i + 1]; k++) {
                    const TIndex j = m_ja[k];
                    m_indices[i][j] = k;
                }
            }
    }
    }

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

            if (it == m_indices[i].end()) {
                return -1;
            }

            return it->second;
        } else {
            const auto lower = m_ja.begin() + m_ia[i];
            const auto upper = m_ja.begin() + m_ia[i + 1];

            const auto it = std::lower_bound(lower, upper, j);

            if (*it != j || it == upper) {
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
        const index size_i = TRowMajor ? rows : cols;
        const index size_j = TRowMajor ? cols : rows;

        assert(length(pattern) == size_i);

        m_rows = static_cast<TIndex>(rows);
        m_cols = static_cast<TIndex>(cols);

        m_ia.resize(size_i + 1);

        m_ia[0] = 0;

        for (TIndex i = 0; i < size_i; i++) {
            const TIndex n = static_cast<TIndex>(pattern[i].size());

            m_ia[i + 1] = m_ia[i] + n;
        }

        m_ja.resize(nb_nonzeros());

        auto ja_it = m_ja.begin();

        for (TIndex i = 0; i < size_i; i++) {
            const TIndex n = static_cast<TIndex>(pattern[i].size());

            for (const auto j : pattern[i]) {
                assert(j < size_j);
                *ja_it++ = static_cast<TIndex>(j);
            }
        }

        if (TIndexMap) {
            m_indices.resize(size_i);

            for (index i = 0; i < size_i; i++) {
                m_indices[i].set_empty_key(-1);
                m_indices[i].resize(m_ia[i + 1] - m_ia[i]);
                for (TIndex k = m_ia[i]; k < m_ia[i + 1]; k++) {
                    const TIndex j = m_ja[k];
                    m_indices[i][j] = k;
                }
            }
        }
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
        
        SparseStructure result;
        result.set(m_rows, m_cols, pattern);

        std::vector<index> value_indices;
        value_indices.reserve(nb_nonzeros() * 2 - n);

        for (const auto& indices_row : indices) {
            value_indices.insert(value_indices.end(), indices_row.begin(), indices_row.end());
        }

        return {result, value_indices};
    }

    /*
    * https://github.com/scipy/scipy/blob/3b36a574dc657d1ca116f6e230be694f3de31afc/scipy/sparse/sparsetools/csr.h#L380
    */
    static Type convert_from(SparseStructure<TScalar, TIndex, !TRowMajor> other, Ref<Vector> values)
    {
        const auto nb_nonzeros = other.nb_nonzeros();

        const auto n = TRowMajor ? other.cols() : other.rows();
        const auto m = TRowMajor ? other.rows() : other.cols();

        std::vector<TIndex> ia(m + 1);
        std::vector<TIndex> ja(nb_nonzeros);

        std::fill(ia.begin(), ia.end(), 0);

        for (TIndex k = 0; k < nb_nonzeros; k++) {
            ia[other.ja(k)] += 1;
        }

        TIndex cumsum = 0;

        for (TIndex j = 0; j < m; j++) {
            const auto temp  = ia[j];
            ia[j] = cumsum;
            cumsum += temp;
        }

        ia[m] = nb_nonzeros;

        Vector a_values = values;

        for (TIndex i = 0; i < n; i++){
            for(TIndex k = other.ia(i); k < other.ia(i + 1); k++){
                const auto j = other.ja(k);
                const auto dest = ia[j];

                ja[dest] = i;
                values[dest] = a_values[k];

                ia[j]++;
            }
        }

        TIndex last = 0;

        for (TIndex j = 0; j <= m; j++){
            const auto temp = ia[j];
            ia[j] = last;
            last = temp;
        }

        return Type(other.rows(), other.cols(), ia, ja);
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

public: // python
    template <typename TModule>
    static void register_python(TModule& m, const std::string& name)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Holder = Pointer<Type>;

        py::class_<Type, Holder>(m, name.c_str())
            // constructors
            .def(py::init<TIndex, TIndex, std::vector<TIndex>, std::vector<TIndex>>(), "rows"_a, "cols"_a, "ia"_a, "ja"_a)
            // static methods
            .def_static("convert_from", &Type::convert_from, "other"_a, "values"_a)
            // methods
            .def("to_general", &Type::to_general)
            // read-only properties
            .def_property_readonly("rows", &Type::rows)
            .def_property_readonly("cols", &Type::cols)
            .def_property_readonly("nb_nonzeros", &Type::nb_nonzeros)
            .def_property_readonly("density", &Type::density)
            .def_property_readonly("ia", py::overload_cast<>(&Type::ia, py::const_))
            .def_property_readonly("ja", py::overload_cast<>(&Type::ja, py::const_));
    }
};

} // namespace eqlib