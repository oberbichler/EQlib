#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <fmt/format.h>

#include <pybind11/pybind11.h>

#include <string>
#include <variant>

namespace EQlib {

template <typename T>
using Pointer = std::shared_ptr<T>;

template <typename T>
using Unique = std::unique_ptr<T>;

template <typename T, typename... TArgs>
Unique<T> new_(TArgs&&... args)
{
    return Unique<T>(new T(std::forward<TArgs>(args)...));
}

using Vector3D = Eigen::Vector3d;

using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Sparse = Eigen::SparseMatrix<double, Eigen::ColMajor>;

template <int N>
using VectorN = Eigen::Matrix<double, N, 1>;

template <typename T>
using Ref = Eigen::Ref<T>;

template <typename T>
using Map = Eigen::Map<T>;

template <typename... Args>
std::string format(Args&&... args)
{
    return fmt::format(std::forward<Args>(args)...);
}

} // namespace