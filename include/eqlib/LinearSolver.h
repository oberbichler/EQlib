#pragma once

#include "Define.h"

#include <string>

namespace eqlib {

class LinearSolver {
private: // types
    using Type = LinearSolver;

private: // variables
    std::string m_solver_name;

public: // constructors
    virtual ~LinearSolver() = default;

public: // methods
    std::string solver_name() const
    {
        return m_solver_name;
    }

    void set_solver_name(const std::string& value)
    {
        m_solver_name = value;
    }

    virtual bool analyze(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a)
    {
        return false;
    }

    virtual bool factorize(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a)
    {
        return false;
    }

    virtual bool solve(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a, Ref<const Vector> b, Ref<Vector> x) = 0;
    
public: // python
    template <typename T>
    class PyLinearSolver : public T {
    public: // constructor
        using T::T;

    public: // methods
        virtual bool analyze(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a) override
        {
            pybind11::gil_scoped_acquire acquire;
            PYBIND11_OVERLOAD(bool, T, analyze, ia, ja, a);
        }

        virtual bool factorize(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a) override
        {
            pybind11::gil_scoped_acquire acquire;
            PYBIND11_OVERLOAD(bool, T, factorize, ia, ja, a);
        }

        virtual bool solve(const std::vector<int>& ia, const std::vector<int>& ja, Ref<const Vector> a, Ref<const Vector> b, Ref<Vector> x) override
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
            // read-only properties
            .def_property("solver_name", &Type::solver_name, &Type::set_solver_name)
            // methods
            .def("analyze", &Type::analyze, "ia"_a, "ja"_a, "a"_a)
            .def("factorize", &Type::factorize, "ia"_a, "ja"_a, "a"_a)
            .def("solve", &Type::solve, "ia"_a, "ja"_a, "a"_a, "b"_a, "x"_a);
    }
};

} // namespace eqlib