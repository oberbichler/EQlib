#pragma once

#include "Define.h"
#include "Parameter.h"
#include "Element.h"
#include "Log.h"
#include "Settings.h"
#include "Timer.h"

#include "LinearSolvers/EigenLinearSolver.h"
#include "LinearSolvers/LinearSolver.h"

#include <sparsehash/dense_hash_map>

#include <tbb/tbb.h>

#include <tsl/hopscotch_set.h>
#include <tsl/robin_set.h>

#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace EQlib {

template <bool TSymmetric = true>
class System
{
private:    // types
    struct Index
    {
        int local;
        int global;

        bool operator<(const Index& other) const noexcept
        {
            return this->global < other.global;
        }
    };

private:    // variables
    std::vector<Pointer<Parameter>> m_dofs;

    google::dense_hash_map<Pointer<Parameter>, int> m_dof_indices;

    int m_nb_free_dofs;
    int m_nb_fixed_dofs;

    int m_max_element_size;

    std::vector<Pointer<Element>> m_elements;
    std::vector<std::vector<Index>> m_index_table;

    std::vector<std::vector<int>> m_pattern;
    Eigen::VectorXi m_col_nonzeros;

    double m_f;
    Vector m_g;
    Sparse m_h;

    Vector m_x;
    Vector m_target;
    Vector m_residual;

    int m_stopping_reason;

    Unique<LinearSolver> m_solver;

    double m_load_factor;

public:     // constructors
    System(
        std::vector<Pointer<Element>> elements,
        Settings linear_solver)
    : m_elements(std::move(elements))
    , m_load_factor(1)
    {
        initialize(linear_solver);
    }

private:    // methods
    void initialize(Settings linear_solver)
    {
        Log::info(1, "==> Initialize system...");
        Log::info(2, "The system consists of {} elements", m_elements.size());

        Timer timer;

        // get dofs

        Log::info(3, "Getting element dofs...");

        const auto nb_elements = m_elements.size();

        std::vector<std::vector<Pointer<Parameter>>> element_dofs(nb_elements);

        m_max_element_size = 0;

        for (size_t i = 0; i < nb_elements; i++) {
            const auto& element = m_elements[i];

            element_dofs[i] = element->dofs();

            const auto element_size = element_dofs[i].size();

            if (element_size > m_max_element_size) {
                m_max_element_size = element_size;
            }
        }

        // create set of unique dofs

        Log::info(3, "Creating set of unique dofs...");

        tsl::robin_set<Pointer<Parameter>> dof_set;
        std::vector<Pointer<Parameter>> free_dofs;
        std::vector<Pointer<Parameter>> fixed_dofs;

        for (size_t i = 0; i < nb_elements; i++) {
            const auto& element = m_elements[i];

            for (const auto& dof : element_dofs[i]) {
                if (dof_set.find(dof) != dof_set.end()) {
                    continue;
                }

                dof_set.insert(dof);

                if (dof->is_fixed()) {
                    fixed_dofs.push_back(dof);
                } else {
                    free_dofs.push_back(dof);
                }
            }
        }

        // store system size

        m_nb_free_dofs = static_cast<int>(free_dofs.size());
        m_nb_fixed_dofs = static_cast<int>(fixed_dofs.size());

        Log::info(2, "The system has {} free and {} fixed dofs", m_nb_free_dofs,
            m_nb_fixed_dofs);

        // create a vector of unique dofs

        m_dofs.reserve(m_nb_free_dofs + m_nb_fixed_dofs);

        m_dofs.insert(m_dofs.end(), free_dofs.begin(), free_dofs.end());
        m_dofs.insert(m_dofs.end(), fixed_dofs.begin(), fixed_dofs.end());

        // create a {dof -> index} map

        Log::info(3, "Creating dof indices...");

        m_dof_indices.resize(m_dofs.size());
        m_dof_indices.set_empty_key(Pointer<Parameter>());

        for (int i = 0; i < m_dofs.size(); i++) {
            m_dof_indices[m_dofs[i]] = i;
        }

        // create index table

        Log::info(3, "Creating index table for elements...");

        for (size_t i = 0; i < nb_elements; i++) {
            const auto& element = m_elements[i];

            const auto& dofs = element_dofs[i];

            const auto nb_dofs = dofs.size();

            std::vector<Index> dof_indices(nb_dofs);

            for (int local = 0; local < nb_dofs; local++) {
                const auto dof = dofs[local];

                int global = free_dofs.size();

                const auto it = m_dof_indices.find(dof);

                if (it != m_dof_indices.end()) {
                    global = it->second;
                }

                dof_indices[local] = {local, global};
            }

            std::sort(dof_indices.begin(), dof_indices.end());

            m_index_table.push_back(dof_indices);
        }

        // analyze pattern

        Log::info(3, "Analyzing pattern...");

        std::vector<tsl::robin_set<int>> pattern(m_nb_free_dofs);

        for (const auto& dof_indices : m_index_table) {
            const size_t nb_dofs = dof_indices.size();

            for (size_t row = 0; row < nb_dofs; row++) {
                const auto row_index = dof_indices[row];

                if (row_index.global >= nb_free_dofs()) {
                    continue;
                }

                for (size_t col = TSymmetric ? row : 0; col < nb_dofs; col++) {
                    const auto col_index = dof_indices[col];

                    if (col_index.global >= nb_free_dofs()) {
                        continue;
                    }

                    pattern[col_index.global].insert(row_index.global);
                }
            }
        }

        m_col_nonzeros = Eigen::VectorXi(m_nb_free_dofs);

        for (int i = 0; i < pattern.size(); i++) {
            m_col_nonzeros[i] = static_cast<int>(pattern[i].size());
        }

        // store data

        m_pattern = std::vector<std::vector<int>>(pattern.size());

        for (int col = 0; col < pattern.size(); col++) {
            std::vector<int> rows;
            rows.insert(rows.end(), pattern[col].begin(), pattern[col].end());
            std::sort(rows.begin(), rows.end());
            m_pattern[col] = std::move(rows);
        }

        // setup system vectors and matrices

        Log::info(3, "Allocating memory...");

        m_h = Sparse(nb_free_dofs(), nb_free_dofs());

        if (nb_free_dofs() > 0) {
            m_h.reserve(m_col_nonzeros);

            for (int col = 0; col < m_pattern.size(); col++) {
                for (const int row : m_pattern[col]) {
                    m_h.insert(row, col);
                }
            }
        }

        Log::info(2, "The system matrix has {} nonzero entries ({:.3f}%)",
            m_h.nonZeros(), m_h.nonZeros() * 100.0 / m_h.size());

        m_g = Vector(nb_free_dofs());

        m_x = Vector(nb_free_dofs());

        m_target = Vector(nb_free_dofs());
        m_residual = Vector(nb_free_dofs());

        // setup solver

        Log::info(3, "Initializing solver...");

        const auto linear_solver_type =
            get_or_default(linear_solver, "type", "pardiso_ldlt");

        if (linear_solver_type == "pardiso_ldlt") {
            Log::info(2, "Using Pardiso LDL^T solver");

            m_solver = new_<PardisoLDLT>();
        } else if (linear_solver_type == "pardiso_llt") {
            Log::info(2, "Using Pardiso LL^T solver");

            m_solver = new_<PardisoLLT>();
        } else if (linear_solver_type == "pardiso_lu") {
            Log::info(2, "Using Pardiso LU solver");

            m_solver = new_<PardisoLU>();
        } else if (linear_solver_type == "conjugate_gradient") {
            Log::info(2, "Using Eigen Conjugate Gradient solver");

            if (TSymmetric) {
                m_solver = new_<SymmetricConjugateGradient>();
            } else {
                m_solver = new_<ConjugateGradient>();
            }
        } else if (linear_solver_type == "sparse_lu") {
            Log::info(2, "Using Eigen Sparse LU solver");

            m_solver = new_<SparseLU>();
        } else {
            throw std::runtime_error("Unknown linear solver");
        }

        Log::info(1, "System initialized in {:.3f} sec", timer.ellapsed());
    }

public:     // getters and setters
    const Pointer<Parameter>& dof(const int index) const
    {
        return m_dofs[index];
    }

    const std::vector<Pointer<Parameter>>& dofs() const
    {
        return m_dofs;
    }

    int nb_dofs() const
    {
        return static_cast<int>(m_dofs.size());
    }

    int nb_elements() const
    {
        return static_cast<int>(m_elements.size());
    }

    int nb_free_dofs() const
    {
        return m_nb_free_dofs;
    }

    int nb_fixed_dofs() const
    {
        return m_nb_fixed_dofs;
    }

    double f() const
    {
        return m_f;
    }

    Vector g() const
    {
        return m_g;
    }

    double g(const int index) const
    {
        return m_g(index);
    }

    const Sparse& h() const
    {
        return m_h;
    }

    double h(const int row, const int col) const
    {
        return m_h.coeff(row, col);
    }

    int h_nb_nonzeros() const
    {
        return static_cast<int>(m_h.nonZeros());
    }

    Vector h_inv_v(Ref<const Vector> v)
    {
        if (nb_free_dofs() > 0) {
            m_solver->factorize(m_h);

            if (m_solver->info() != Eigen::Success) {
                throw std::runtime_error("Factorization failed");
            }

            Vector x(nb_free_dofs());

            m_solver->solve(v, x);

            if (m_solver->info() != Eigen::Success) {
                throw std::runtime_error("Solve failed");
            }

            return x;
        } else {
            return Vector(0);
        }
    }

    Vector h_v(Ref<const Vector> v) const
    {
        if (TSymmetric) {
        return m_h.selfadjointView<Eigen::Upper>() * v;
        } else {
            return m_h * v;
    }
    }

    Vector delta() const
    {
        Vector result(nb_free_dofs());

        for (index i = 0; i < length(result); i++) {
            result(i) = dof(i)->delta();
        }

        return result;
    }

    void set_delta(Ref<const Vector> value) const
    {
        if (length(value) != nb_free_dofs()) {
            throw std::runtime_error("Invalid size");
        }

        for (index i = 0; i < length(value); i++) {
            dof(i)->set_delta(value[i]);
        }
    }

    void set_delta(double* const value) const
    {
        set_delta(Map<const Vector>(value, nb_free_dofs()));
    }

    Vector x() const
    {
        Vector result(nb_free_dofs());

        for (int i = 0; i < result.size(); i++) {
            result[i] = m_dofs[i]->act_value();
        }

        return result;
    }

    void set_x(Ref<const Vector> value) const
    {
        if (value.size() != nb_free_dofs()) {
            throw std::runtime_error("Invalid size");
        }

        for (int i = 0; i < value.size(); i++) {
            m_dofs[i]->set_act_value(value[i]);
        }
    }

    void set_x(double* const values) const
    {
        set_x(Map<const Vector>(values, nb_free_dofs()));
    }

    Vector residual() const
    {
        return m_residual;
    }

    double load_factor() const
    {
        return m_load_factor;
    }

    void set_load_factor(const double value)
    {
        m_load_factor = value;
    }

    std::vector<Pointer<Element>> elements() const
    {
        return m_elements;
    }

    std::vector<int> element_indices(int index) const
    {
        const auto& indices = m_index_table[index];

        std::vector<int> element_indices(indices.size());

        for (const auto& index : indices) {
            element_indices[index.local] = index.global;
        }

        return element_indices;
    }

    std::string message() const
    {
        switch (m_stopping_reason) {
        case -1:
            return "Not solved";
        case 0:
            return "A solution was found, given rtol";
        case 1:
            return "A solution was found, given xtol";
        case 2:
            return "The iteration limit was reached";
        default:
            return "Error. Unknown stopping reason";
        }
    }

public:     // methods
    void add_diagonal(const double value)
    {
        for (int i = 0; i < m_h.rows(); i++) {
            m_h.coeffRef(i, i) += value;
        }
    }

    int dof_index(const Pointer<Parameter>& dof) const
    {
        const auto it = m_dof_indices.find(dof);

        if (it == m_dof_indices.end()) {
            return -1;
        }

        return it->second;
    }

    template<int TOrder>
    void assemble(const bool parallel)
    {
        static_assert(0 <= TOrder && TOrder <= 2);

        assemble<TOrder>(parallel, m_f, m_g, m_h);
    }

    template<int TOrder>
    void assemble(const bool parallel, double& f, Ref<Vector> g, Ref<Sparse> h)
    {
        static_assert(0 <= TOrder && TOrder <= 2);

        auto begin = tbb::make_zip_iterator(m_elements.begin(),
            m_index_table.begin());
        auto end = tbb::make_zip_iterator(m_elements.end(),
            m_index_table.end());

        if (parallel) {
            switch (TOrder) {
            case 0:
                assemble_parallel<0>(begin, end, f, g, h);
                break;
            case 1:
                assemble_parallel<1>(begin, end, f, g, h);
                break;
            case 2:
                assemble_parallel<2>(begin, end, f, g, h);
                break;
            }
        } else {
            switch (TOrder) {
            case 0:
                assemble_serial<0>(begin, end, f, g, h);
                break;
            case 1:
                assemble_serial<1>(begin, end, f, g, h);
                break;
            case 2:
                assemble_serial<2>(begin, end, f, g, h);
                break;
            }
        }
    }

    template <int TOrder, typename TIterator>
    void assemble_parallel(TIterator begin, TIterator end, double& f,
        Ref<Vector> g, Ref<Sparse> h, const bool init_zero = true)
    {
        static_assert(0 <= TOrder && TOrder <= 2);

        pybind11::gil_scoped_release release;

        // compute g and h

        tbb::combinable<double> c_f(0);
        tbb::combinable<Vector> c_g(Vector::Zero(m_g.size()));
        tbb::combinable<Vector> c_h(Vector::Zero(m_h.nonZeros()));

        tbb::combinable<Vector> buffer_g([=]() {
            return Vector(m_max_element_size);
        });
        tbb::combinable<Matrix> buffer_h([=]() {
            return Matrix(m_max_element_size, m_max_element_size);
        });

        Vector dummy_vector = Vector(0);
        Matrix dummy_matrix = Matrix(0, 0);
        Sparse dummy_sparse = Sparse(0, 0);

        tbb::parallel_for(tbb::blocked_range<decltype(begin)>(begin, end, 128),
            [&](const tbb::blocked_range<decltype(begin)> &range) {
            Log::info(10, "New task with {} items", range.size());

            double& local_f = c_f.local();

            if (TOrder == 0) {
                auto& local_g = dummy_vector;
                auto& local_h = dummy_sparse;

                auto& local_buffer_g = dummy_vector;
                auto& local_buffer_h = dummy_matrix;

                assemble_serial<TOrder>(range.begin(), range.end(),
                    local_buffer_g, local_buffer_h, local_f, local_g, local_h,
                    false);
            } else if (TOrder == 1) {
                auto& local_g = c_g.local();
                auto& local_h = dummy_sparse;

                auto& local_buffer_g = buffer_g.local();
                auto& local_buffer_h = dummy_matrix;

                assemble_serial<TOrder>(range.begin(), range.end(),
                    local_buffer_g, local_buffer_h, local_f, local_g, local_h,
                    false);
            } else if (TOrder == 2) {
                auto& local_g = c_g.local();
                auto local_h = Map<Sparse>(m_h.rows(), m_h.cols(),
                    m_h.nonZeros(), m_h.outerIndexPtr(), m_h.innerIndexPtr(),
                    c_h.local().data());

                auto& local_buffer_g = buffer_g.local();
                auto& local_buffer_h = buffer_h.local();

                assemble_serial<TOrder>(range.begin(), range.end(),
                    local_buffer_g, local_buffer_h, local_f, local_g, local_h,
                    false);
            }
        });

        const auto sum_f = c_f.combine(std::plus<double>());

        if (init_zero) {
            f = sum_f;
        } else {
            f += sum_f;
        }

        if (TOrder > 0) {
            const auto sum_g = c_g.combine(std::plus<Vector>());

            if (init_zero) {
                g = sum_g;
            } else {
                g += sum_g;
            }
        }

        if (TOrder > 1) {
            const auto sum_h = c_h.combine(std::plus<Vector>());

            if (init_zero) {
                Map<Vector>(h.valuePtr(), h.nonZeros()) = sum_h;
            } else {
                Map<Vector>(h.valuePtr(), h.nonZeros()) += sum_h;
            }
        }
    }

    template <int TOrder, typename TIterator>
    void assemble_serial(TIterator begin, TIterator end, double& f,
        Ref<Vector> g, Ref<Sparse> h, const bool init_zero = true)
    {
        static_assert(0 <= TOrder && TOrder <= 2);

        Vector buffer_g(m_max_element_size);
        Matrix buffer_h(m_max_element_size, m_max_element_size);

        assemble_serial<TOrder>(begin, end, buffer_g, buffer_h, f, g, h,
            init_zero);
    }

    template <int TOrder, typename TIterator>
    void assemble_serial(TIterator begin, TIterator end, Ref<Vector> buffer_g,
        Ref<Matrix> buffer_h, double& f, Ref<Vector> g, Ref<Sparse> h,
        const bool init_zero = true)
    {
        static_assert(0 <= TOrder && TOrder <= 2);

        if (init_zero) {
            f = 0;

            if (TOrder > 0) {
                g.setZero();
            }

            if (TOrder > 1) {
                Map<Vector>(h.valuePtr(), h.nonZeros()).setZero();
            }
        }

        for (auto it = begin; it != end; ++it) {
            const auto& element = std::get<0>(*it);
            const auto& dof_indices = std::get<1>(*it);

            const size_t nb_dofs = dof_indices.size();

            Ref<Vector> element_g = buffer_g.head(TOrder > 0 ? nb_dofs : 0);
            Ref<Matrix> element_h = buffer_h.topLeftCorner(TOrder > 1 ?
                nb_dofs : 0, TOrder > 1 ? nb_dofs : 0);

            const auto element_f = element->compute(element_g, element_h);

            f += element_f;

            if (TOrder < 1) {
                continue;
            }

            for (size_t row = 0; row < nb_dofs; row++) {
                const auto row_index = dof_indices[row];

                if (row_index.global >= g.size()) {
                    continue;
                }

                g(row_index.global) += element_g(row_index.local);

                if (TOrder < 2) {
                    continue;
                }

                for (size_t col = TSymmetric ? row : 0; col < nb_dofs; col++) {
                    const auto col_index = dof_indices[col];

                    if (col_index.global >= g.size()) {
                        continue;
                    }

                    h.coeffRef(row_index.global, col_index.global) +=
                        element_h(row_index.local, col_index.local);
                }
            }
        }
    }

    void compute(const int order, const bool parallel)
    {
        Log::info(1, "==> Computing system...");

        Timer timer;

        switch (order) {
        case 0:
            assemble<0>(parallel, m_f, m_g, m_h);
            break;
        case 1:
            assemble<1>(parallel, m_f, m_g, m_h);
            break;
        case 2:
            assemble<2>(parallel, m_f, m_g, m_h);
            break;
        default:
            throw std::runtime_error("Invalid order");
        }

        Log::info(1, "System computed in {:.3f} sec", timer.ellapsed());
    }

    void solve(const int maxiter, const double rtol, const double xtol,
        const double damping, const bool parallel)
    {
        // setup

        Log::info(1, "==> Solving nonlinear system...");

        Timer timer;

        for (int i = 0; i < nb_free_dofs(); i++) {
            m_target[i] = m_dofs[i]->target();
        }

        if (m_load_factor != 1.0) {
            m_target *= m_load_factor;
        }

        for (int iteration = 0; ; iteration++) {
            // check max iterations

            if (iteration >= maxiter) {
                m_stopping_reason = 2;
                Log::info(2, "Stopped because iteration >= {}", maxiter);
                break;
            }

            Log::info(2, "Iteration {}", iteration + 1);

            // compute g and h

            Log::info(2, "Computing system...");

            assemble<2>(parallel, m_f, m_g, m_h);

            Log::info(2, "The current value is {}", m_f);

            // check residual

            Log::info(2, "Computing residual...");

            m_residual = m_target - m_g;

            const double rnorm = m_residual.norm();

            Log::info(2, "The norm of the residual is {}", rnorm);

            // check residual norm

            if (rnorm < rtol) {
                m_stopping_reason = 0;
                Log::info(2, "Stopped because rnorm < {}", rtol);
                break;
            }

            // solve iteration

            Log::info(2, "Solving the linear equation system...");

            if (damping != 0.0) {
                add_diagonal(damping);
            }

            m_x = h_inv_v(m_residual);

            // update system

            Log::info(2, "Updating system...");

            for (int i = 0; i < nb_free_dofs(); i++) {
                m_dofs[i]->set_delta(m_dofs[i]->delta() + m_x(i));
            }

            // check x norm

            const double xnorm = m_x.norm();

            Log::info(2, "The norm of the step is {}", xnorm);

            if (xnorm < xtol) {
                Log::info(2, "Stopped because xnorm < {}", xtol);
                m_stopping_reason = 1;
                break;
            }
        }

        for (int i = 0; i < nb_free_dofs(); i++) {
            m_dofs[i]->set_residual(m_residual(i));
        }

        switch (m_stopping_reason) {
        case 2:
            Log::warn("The maximum number of iterations has been reached");
            break;
        case 3:
            Log::error("An unknown error has occurred");
            break;
        }

        Log::info(1, "System solved in {:.3f} sec", timer.ellapsed());
    }

    void solve_linear(const bool parallel, const bool update_dofs)
    {
        // setup

        Log::info(1, "==> Solving linear system...");

        Timer timer;

        for (int i = 0; i < nb_free_dofs(); i++) {
            m_target[i] = m_dofs[i]->target();
        }

        if (m_load_factor != 1.0) {
            m_target *= m_load_factor;
        }

        m_residual = m_target - m_g;

        // compute g and h

        assemble<2>(parallel, m_f, m_g, m_h);

        // solve

        m_x = h_inv_v(m_residual);

        // update system

        if (update_dofs) {
            for (int i = 0; i < nb_free_dofs(); i++) {
                m_dofs[i]->set_delta(m_x(i));
            }
        }

        Log::info(1, "System solved in {:.3f} sec", timer.ellapsed());
    }

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Type = System<TSymmetric>;
        using Holder = Pointer<Type>;

        const std::string name = TSymmetric ? "SymmetricSystem" : "System";

        py::class_<Type, Holder>(m, name.c_str())
            // constructors
            .def(py::init<std::vector<Pointer<EQlib::Element>>, Settings>(),
                "elements"_a, "linear_solver"_a = Settings())
            // properties
            .def_property("load_factor", &Type::load_factor,
                &Type::set_load_factor)
            .def_property("delta", &Type::delta,
                py::overload_cast<Ref<const Vector>>(&Type::set_delta, py::const_))
            .def_property("x", &Type::x,
                py::overload_cast<Ref<const Vector>>(&Type::set_x, py::const_))
            // readonly properties
            .def_property_readonly("dofs", &Type::dofs)
            .def_property_readonly("nb_dofs", &Type::nb_dofs)
            .def_property_readonly("nb_elements", &Type::nb_elements)
            .def_property_readonly("elements", &Type::elements)
            .def_property_readonly("f", &Type::f)
            .def_property_readonly("g",
                py::overload_cast<void>(&Type::g, py::const_))
            .def_property_readonly("h",
                py::overload_cast<>(&Type::h, py::const_))
            .def_property_readonly("message", &Type::message)
            .def_property_readonly("nb_free_dofs", &Type::nb_free_dofs)
            .def_property_readonly("nb_fixed_dofs", &Type::nb_fixed_dofs)
            .def_property_readonly("residual", &Type::residual)
            // methods
            .def("add_diagonal", &Type::add_diagonal, "value"_a)
            .def("compute", &Type::compute, "order"_a=2, "parallel"_a=true)
            .def("dof_index", &Type::dof_index)
            .def("element_indices", &Type::element_indices, "index"_a)
            .def("h_inv_v", &Type::h_inv_v)
            .def("h_v", &Type::h_v)
            .def("solve", &Type::solve, "maxiter"_a = 100, "rtol"_a = 1e-7,
                "xtol"_a = 1e-7, "damping"_a=0.0, "parallel"_a=true)
            .def("solve_linear", &Type::solve_linear, "parallel"_a=true,
                "update_dofs"_a=true)
        ;
    }
};

} // namespace EQlib