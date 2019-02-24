#pragma once

#include <Eigen/Core>
#include <Eigen/PardisoSupport>
#include <Eigen/Sparse>

#include <pybind11/pybind11.h>

using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Sparse = Eigen::SparseMatrix<double, Eigen::ColMajor>;

#if defined EIGEN_USE_MKL_ALL
using SparseSolver = Eigen::PardisoLLT<Sparse, Eigen::Upper>;
#else
using SparseSolver = Eigen::SparseLU<Sparse>;
#endif

namespace py = pybind11;

namespace EQlib {

template <typename T>
T get_or_default(py::dict options, std::string key, T default_value) {
    if (!options.contains(key.c_str())) {
        return default_value;
    }
    return options[key.c_str()].cast<T>();
}

}