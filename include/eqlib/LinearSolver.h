#pragma once

#include "Define.h"

#include <string>

namespace eqlib {

class LinearSolver {
private: // types
    using Type = LinearSolver;

public: // constructors
    virtual ~LinearSolver() = default;

public: // methods
    virtual std::string solver_name() const = 0;

    virtual bool analyze(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a) = 0;

    virtual bool factorize(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a) = 0;

    virtual bool solve(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a, Ref<const Vector> b, Ref<Vector> x) = 0;
    
public: // python
    template <typename T>
    class PyLinearSolver : public T {
    public: // constructor
        using T::T;

    public: // methods
        virtual std::string solver_name() const
        {
            pybind11::gil_scoped_acquire acquire;
            PYBIND11_OVERLOAD_PURE(std::string, T, solver_name);
        }
    
        virtual bool analyze(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a)
        {
            pybind11::gil_scoped_acquire acquire;
            PYBIND11_OVERLOAD_PURE(bool, T, analyze, ia, ja, a);
        }

        virtual bool factorize(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a)
        {
            pybind11::gil_scoped_acquire acquire;
            PYBIND11_OVERLOAD_PURE(bool, T, factorize, ia, ja, a);
        }

        virtual bool solve(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a, Ref<const Vector> b, Ref<Vector> x)
        {
            pybind11::gil_scoped_acquire acquire;
            PYBIND11_OVERLOAD_PURE(bool, T, factorize, ia, ja, a, b, x);
        }
    };

    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Trampoline = PyLinearSolver<Type>;
        using Holder = Pointer<Type>;

        py::class_<Type, Trampoline, Holder>(m, "LinearSolver")
            .def(py::init<>())
            // methods
            .def("solver_name", &Type::solver_name)
            .def("analyze", &Type::analyze, "ia"_a, "ja"_a, "a"_a)
            .def("factorize", &Type::factorize, "ia"_a, "ja"_a, "a"_a)
            .def("solve", &Type::solve, "ia"_a, "ja"_a, "a"_a, "b"_a, "x"_a);
    }
};

} // namespace eqlib