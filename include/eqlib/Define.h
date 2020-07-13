#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>

#include <fmt/format.h>

#include <pybind11/pybind11.h>

#include <tsl/robin_map.h>
#include <tsl/robin_set.h>

#include <sparsehash/dense_hash_map>

#include <limits>
#include <string>
#include <variant>

namespace eqlib {

#if defined(_MSC_VER)
#define EQLIB_INLINE __forceinline
#else
#define EQLIB_INLINE __attribute__((always_inline)) inline
#endif

using index = std::ptrdiff_t;

template <typename T>
EQLIB_INLINE index length(const T& container) noexcept
{
    return static_cast<index>(container.size());
}

// --- memory

template <typename T>
using Pointer = std::shared_ptr<T>;

template <typename T>
using Unique = std::unique_ptr<T>;

template <typename T, typename... TArgs>
Unique<T> new_(TArgs&&... args)
{
    return Unique<T>(new T(std::forward<TArgs>(args)...));
}

// --- hashset

template <typename TKey>
using RobinSet = tsl::robin_set<TKey>;

template <typename TKey, typename TValue>
using RobinMap = tsl::robin_map<TKey, TValue>;

template <typename TKey, typename TValue>
using DenseMap = google::dense_hash_map<TKey, TValue>;

// --- linear algebra

using Vector3D = Eigen::Vector3d;

using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using Vector = Eigen::Matrix<double, 1, Eigen::Dynamic>;
using Sparse = Eigen::SparseMatrix<double, Eigen::RowMajor>;

template <int N>
using VectorN = Eigen::Matrix<double, N, 1>;

template <typename T>
using Ref = Eigen::Ref<T>;

template <typename T>
using Map = Eigen::Map<T>;

// --- format

template <typename... Args>
std::string format(Args&&... args)
{
    return fmt::format(std::forward<Args>(args)...);
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v)
{
    if (!v.empty()) {
        out << '[';
        std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
        out << "\b\b]";
    }
    return out;
}

const double infinity = std::numeric_limits<double>::infinity();

} // namespace eqlib