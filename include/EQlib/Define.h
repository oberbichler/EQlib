#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <fmt/format.h>

#include <pybind11/pybind11.h>

#include <string>

namespace EQlib {

using Vector3D = Eigen::Vector3d;

using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Sparse = Eigen::SparseMatrix<double, Eigen::ColMajor>;

template <typename T>
using Ref = Eigen::Ref<T>;

template <typename T>
using Map = Eigen::Map<T>;

namespace py = pybind11;

template <typename T>
T get_or_default(py::dict options, std::string key, T default_value) {
    if (!options.contains(key.c_str())) {
        return default_value;
    }
    return options[key.c_str()].cast<T>();
}

template <typename... Args>
std::string format(Args&&... args)
{
    return fmt::format(std::forward<Args>(args)...);
}

std::string format_number(const double value)
{
    if (value == std::numeric_limits<double>::min()) {
        return "-inf";
    }

    if (value == std::numeric_limits<double>::max()) {
        return "inf";
    }

    return format("{}", value);
}

} // namespace EQlib