#pragma once

#include "Define.h"
#include "Constraint.h"
#include "LinearSolver.h"
#include "Objective.h"
#include "ProblemData.h"
#include "Settings.h"
#include "Timer.h"
#include "SparseStructure.h"

#include <sparsehash/dense_hash_map>

#include <tbb/tbb.h>

#include <tsl/robin_set.h>

#include <set>
#include <tuple>
#include <utility>
#include <vector>

namespace EQlib {

class Problem
{
private:    // types
    using ElementsF = std::vector<Pointer<Objective>>;
    using ElementsG = std::vector<Pointer<Constraint>>;

    using Equations = std::vector<Pointer<Equation>>;
    using Variables = std::vector<Pointer<Variable>>;

    struct Index
    {
        index local;
        index global;

        Index(index local, index global)
        : local(local), global(global)
        {
        }

        bool operator<(const Index& other) const noexcept
        {
            return global < other.global;
        }
    };

private:    // variables
    double m_sigma;

    bool m_parallel;

    ElementsF m_elements_f;
    ElementsG m_elements_g;

    std::vector<Pointer<Equation>> m_equations;
    std::vector<Pointer<Variable>> m_variables;

    google::dense_hash_map<Pointer<Equation>, index> m_equation_indices;
    google::dense_hash_map<Pointer<Variable>, index> m_variable_indices;

    std::vector<index> m_element_f_nb_variables;
    std::vector<index> m_element_g_nb_variables;
    std::vector<index> m_element_g_nb_equations;

    index m_max_element_n;
    index m_max_element_m;

    std::vector<std::vector<Index>> m_element_f_variable_indices;
    std::vector<std::vector<Index>> m_element_g_equation_indices;
    std::vector<std::vector<Index>> m_element_g_variable_indices;

    SparseStructure<double, int, false> m_dg_structure;
    SparseStructure<double, int, false> m_hl_structure;

    ProblemData m_data;
    tbb::combinable<ProblemData> m_local_data;

    LinearSolver m_linear_solver;

public:     // constructors
    Problem(
        ElementsF elements_f,
        ElementsG elements_g,
        Settings linear_solver)
    : m_elements_f(std::move(elements_f))
    , m_elements_g(std::move(elements_g))
    , m_sigma(1.0)
    , m_parallel(false)
    , m_max_element_n(0)
    , m_max_element_m(0)
    {
        Log::info(1, "==> Initialize problem...");

        Timer timer;

        const auto nb_elements_f = length(m_elements_f);
        const auto nb_elements_g = length(m_elements_g);

        Log::info(2, "The objective consists of {} elements", nb_elements_f);
        Log::info(2, "The constraints consist of {} elements", nb_elements_g);


        Log::info(3, "Getting equations and variables...");

        m_element_f_nb_variables.resize(nb_elements_f);
        m_element_g_nb_variables.resize(nb_elements_g);
        m_element_g_nb_equations.resize(nb_elements_g);

        std::vector<Variables> variables_f(nb_elements_f);

        std::vector<Equations> equations_g(nb_elements_g);
        std::vector<Variables> variables_g(nb_elements_g);

        for (index i = 0; i < nb_elements_f; i++) {
            const auto& element = *m_elements_f[i];

            variables_f[i] = element.variables();

            const index nb_variables = length(variables_f[i]);

            m_element_f_nb_variables[i] = nb_variables;
            m_max_element_n = std::max(m_max_element_n, nb_variables);
        }

        for (index i = 0; i < nb_elements_g; i++) {
            const auto& element = *m_elements_g[i];

            equations_g[i] = element.equations();
            variables_g[i] = element.variables();

            const index nb_equations = length(equations_g[i]);
            const index nb_variables = length(variables_g[i]);

            m_element_g_nb_variables[i] = nb_variables;
            m_element_g_nb_equations[i] = nb_equations;

            m_max_element_n = std::max(m_max_element_n, nb_variables);
            m_max_element_m = std::max(m_max_element_m, nb_equations);
        }


        Log::info(3, "Creating the set of unique equations...");

        tsl::robin_set<Pointer<Equation>> equation_set;

        for (const auto& equations : equations_g) {
            for (const auto& equation : equations) {
                if (!equation->is_active()) {
                    continue;
                }

                const auto [_, is_new] = equation_set.insert(equation);

                if (is_new) {
                    m_equations.push_back(equation);
                }
            }
        }


        Log::info(3, "Creating the set of unique variables...");

        tsl::robin_set<Pointer<Variable>> variable_set;

        for (const auto& variables : variables_f) {
            for (const auto& variable : variables) {
                if (!variable->is_active()) {
                    continue;
                }

                const auto [_, is_new] = variable_set.insert(variable);

                if (is_new) {
                    m_variables.push_back(variable);
                }
            }
        }

        for (const auto& variables : variables_g) {
            for (const auto& variable : variables) {
                if (!variable->is_active()) {
                    continue;
                }

                const auto [_, is_new] = variable_set.insert(variable);

                if (is_new) {
                    m_variables.push_back(variable);
                }
            }
        }

        const auto nb_equations = length(m_equations);
        const auto nb_variables = length(m_variables);

        Log::info(2, "The problem contains {} variables", nb_variables);
        Log::info(2, "The problem contains {} constraint equations", nb_equations);


        Log::info(3, "Compute indices for variables and equations...");

        m_equation_indices.set_empty_key(nullptr);
        m_variable_indices.set_empty_key(nullptr);

        m_equation_indices.resize(nb_equations);
        m_variable_indices.resize(nb_variables);

        for (index i = 0; i < length(m_equations); i++) {
            const auto& equation = m_equations[i];
            m_equation_indices[equation] = i;
        }

        for (index i = 0; i < length(m_variables); i++) {
            const auto& variable = m_variables[i];
            m_variable_indices[variable] = i;
        }


        Log::info(3, "Compute indices for elements...");

        // variable indices f

        m_element_f_variable_indices.resize(nb_elements_f);

        for (index i = 0; i < nb_elements_f; i++) {
            const auto& variables = variables_f[i];

            std::vector<Index> variable_indices;
            variable_indices.reserve(variables.size());

            for (index local = 0; local < length(variables); local++) {
                const auto& variable = variables[local];

                if (!variable->is_active()) {
                    continue;
                }

                const auto global = m_variable_indices[variable];

                variable_indices.emplace_back(local, global);
            }

            std::sort(variable_indices.begin(), variable_indices.end());

            m_element_f_variable_indices[i] = std::move(variable_indices);
        }

        // equation indices g

        m_element_g_equation_indices.resize(nb_elements_g);

        for (index i = 0; i < nb_elements_g; i++) {
            const auto& equations = equations_g[i];

            std::vector<Index> equation_indices;
            equation_indices.reserve(equations.size());

            for (index local = 0; local < length(equations); local++) {
                const auto& equation = equations[local];

                if (!equation->is_active()) {
                    continue;
                }

                const auto global = m_equation_indices[equation];

                equation_indices.emplace_back(local, global);
            }

            m_element_g_equation_indices[i] = std::move(equation_indices);
        }

        // variable indices g

        m_element_g_variable_indices.resize(nb_elements_g);

        for (index i = 0; i < nb_elements_g; i++) {
            const auto& variables = variables_g[i];

            std::vector<Index> variable_indices;
            variable_indices.reserve(variables.size());

            for (index local = 0; local < length(variables); local++) {
                const auto& variable = variables[local];

                if (!variable->is_active()) {
                    continue;
                }

                const auto global = m_variable_indices[variable];

                variable_indices.emplace_back(local, global);
            }

            std::sort(variable_indices.begin(), variable_indices.end());

            m_element_g_variable_indices[i] = std::move(variable_indices);
        }


        Log::info(3, "Analyse sparse patterns...");

        const auto n = length(m_variables);
        const auto m = length(m_equations);

        std::vector<std::set<index>> m_pattern_dg(n);
        std::vector<std::set<index>> m_pattern_hl(n);

        for (index i = 0; i < length(m_elements_f); i++) {
            const auto& variable_indices = m_element_f_variable_indices[i];

            for (index col_i = 0; col_i < length(variable_indices); col_i++) {
                const auto col = variable_indices[col_i];

                for (index row_i = col_i; row_i < length(variable_indices); row_i++) {
                    const auto row = variable_indices[row_i];

                    m_pattern_hl[col.global].insert(row.global);
                }
            }
        }

        for (index i = 0; i < length(m_elements_g); i++) {
            const auto& equation_indices = m_element_g_equation_indices[i];
            const auto& variable_indices = m_element_g_variable_indices[i];

            for (const auto row : equation_indices) {
                for (const auto col : variable_indices) {
                    m_pattern_dg[col.global].insert(row.global);
                }
            }

            for (index col_i = 0; col_i < length(variable_indices); col_i++) {
                const auto col = variable_indices[col_i];

                for (index row_i = col_i; row_i < length(variable_indices); row_i++) {
                    const auto row = variable_indices[row_i];

                    m_pattern_hl[col.global].insert(row.global);
                }
            }
        }


        Log::info(3, "Allocate memory...");

        m_dg_structure.set(m, n, m_pattern_dg);
        m_hl_structure.set(n, n, m_pattern_hl);

        m_data.resize(n, m, m_dg_structure.nb_nonzeros(), m_hl_structure.nb_nonzeros(), m_max_element_n, m_max_element_m);


        Log::info(2, "Problem initialized in {} sec", timer.ellapsed());
    }

private:    // methods: computation
    template <index TOrder, typename TBegin, typename TEnd>
    void compute_elements_f(ProblemData& data, const TBegin& begin, const TEnd& end)
    {
        static_assert(0 <= TOrder && TOrder <= 2);

        for (index i = begin; i != end; ++i) {
            const auto& element_f = *m_elements_f[i];

            if (!element_f.is_active()) {
                continue;
            }

            const auto& variable_indices = m_element_f_variable_indices[i];

            const auto n = m_element_f_nb_variables[i];

            Timer timer_element_allocate;

            Map<Vector> g(data.m_buffer.data(), n);
            Map<Matrix> h(data.m_buffer.data() + n, n, n);

            data.m_timer_allocate += timer_element_allocate.ellapsed();

            Timer timer_element_compute;

            const double f = element_f.compute(g, h);

            data.m_timer_compute += timer_element_compute.ellapsed();

            Timer timer_element_assemble;

            data.f() += f;

            if constexpr(TOrder < 1) {
                continue;
            }

            for (index col_i = 0; col_i != length(variable_indices); ++col_i) {
                const auto col = variable_indices[col_i];

                data.df(col.global) += g(col.local);

                if constexpr(TOrder < 2) {
                    continue;
                }

                for (index row_i = col_i; row_i < length(variable_indices); row_i++) {
                    const auto row = variable_indices[row_i];

                    auto hl_values = data.hl();

                    m_hl_structure.coeff_ref(hl_values, row.global, col.global) += h(row.local, col.local);
                }
            }

            data.m_timer_assemble += timer_element_assemble.ellapsed();
        }
    }

    template <index TOrder, typename TBegin, typename TEnd>
    void compute_elements_g(ProblemData& data, const TBegin& begin, const TEnd& end)
    {
        static_assert(0 <= TOrder && TOrder <= 2);

        for (index i = begin; i != end; ++i) {
            const auto& element_g = *m_elements_g[i];

            if (!element_g.is_active()) {
                continue;
            }

            const auto& equation_indices = m_element_g_equation_indices[i];
            const auto& variable_indices = m_element_g_variable_indices[i];

            const auto m = m_element_g_nb_equations[i];
            const auto n = m_element_g_nb_variables[i];

            if (n == 0 || m == 0) { // FIXME: check reduces counts
                continue;
            }

            Timer timer_element_allocate;

            Vector fs(m);
            std::vector<Ref<Vector>> gs;
            std::vector<Ref<Matrix>> hs;

            gs.reserve(m);
            hs.reserve(m);

            for (index k = 0; k < m; k++) {
                Map<Vector> g(data.m_buffer.data() + k * n, n);
                Map<Matrix> h(data.m_buffer.data() + m * n + k * n * n, n, n);
                gs.push_back(g);
                hs.push_back(h);
            }

            data.m_timer_allocate += timer_element_allocate.ellapsed();

            Timer timer_element_compute;

            element_g.compute(fs, gs, hs);

            data.m_timer_compute += timer_element_compute.ellapsed();

            Timer timer_element_assemble;

            for (const auto& equation_index : equation_indices) {
                const auto& equation = m_equations[equation_index.global];

                data.g(equation_index.global) += fs(equation_index.local);

                if constexpr(TOrder < 1) {
                    continue;
                }

                auto& local_g = gs[equation_index.local];
                auto& local_h = hs[equation_index.local];

                local_h *= equation->multiplier();

                for (index col_i = 0; col_i < length(variable_indices); col_i++) {
                    const auto col = variable_indices[col_i];

                    auto dg_values = data.dg();

                    m_dg_structure.coeff_ref(dg_values, equation_index.global, col.global) += local_g(col.local);

                    if constexpr(TOrder < 2) {
                        continue;
                    }

                    for (index row_i = col_i; row_i < length(variable_indices); row_i++) {
                        const auto row = variable_indices[row_i];

                        auto hl_values = data.hl();

                        m_hl_structure.coeff_ref(hl_values, row.global, col.global) += local_h(row.local, col.local);
                    }
                }
            }

            data.m_timer_assemble += timer_element_assemble.ellapsed();
        }
    }

    template <typename TBegin, typename TEnd>
    void compute_elements_f(const index order, ProblemData& data, const TBegin& begin, const TEnd& end)
    {
        switch (order) {
        case 0:
            compute_elements_f<0>(m_data, begin, end);
            break;
        case 1:
            compute_elements_f<1>(m_data, begin, end);
            break;
        case 2:
            compute_elements_f<2>(m_data, begin, end);
            break;
        }
    }

    template <typename TBegin, typename TEnd>
    void compute_elements_g(const index order, ProblemData& data, const TBegin& begin, const TEnd& end)
    {
        switch (order) {
        case 0:
            compute_elements_g<0>(m_data, begin, end);
            break;
        case 1:
            compute_elements_g<1>(m_data, begin, end);
            break;
        case 2:
            compute_elements_g<2>(m_data, begin, end);
            break;
        }
    }

public:     // methods: computation
    void compute(const index order = 2)
    {
        if (order < 0 || 2 < order) {
            throw std::invalid_argument("order");
        }

        Log::info(1, "==> Compute problem...");

        Timer timer;

        m_data.set_zero();

        Log::info(2, "Compute objective...");

        if (!m_parallel) {
            compute_elements_f(order, m_data, 0, nb_elements_f());
        } else {
            tbb::parallel_for(tbb::blocked_range<index>(0, nb_elements_f(), 1000),
                [&](const tbb::blocked_range<index>& range) {
                    Log::info(5, "Launch kernel with {} elements", range.size());
                    auto& local_data = m_local_data.local();
                    local_data.resize(nb_variables(), nb_equations(), m_dg_structure.nb_nonzeros(), m_hl_structure.nb_nonzeros(), m_max_element_n, m_max_element_m);
                    compute_elements_f(order, local_data, range.begin(), range.end());
                }, tbb::static_partitioner()
            );
        }

        m_data.f() *= sigma();

        if (order > 0) {
            m_data.df() *= sigma();
        }

        if (order > 1) {
            m_data.hl() *= sigma();
        }

        Log::info(2, "Compute constraints...");

        if (!m_parallel) {
            compute_elements_g(order, m_data, 0, nb_elements_g());
        } else {
            tbb::parallel_for(tbb::blocked_range<index>(0, nb_elements_g(), 1000),
                [&](const tbb::blocked_range<index>& range) {
                    Log::info(5, "Launch kernel with {} elements", range.size());
                    auto& local_data = m_local_data.local();
                    local_data.resize(nb_variables(), nb_equations(), m_dg_structure.nb_nonzeros(), m_hl_structure.nb_nonzeros(), m_max_element_n, m_max_element_m);
                    compute_elements_g(order, local_data, range.begin(), range.end());
                }, tbb::static_partitioner()
            );

            Log::info(5, "Combine results...");

            m_local_data.combine_each([&](const ProblemData& local) {
                m_data += local;
            });
        }

        Log::info(2, "Problem computed in {} sec", timer.ellapsed());

        Log::info(3, "Elements allocated in {} sec", m_data.m_timer_allocate);
        Log::info(3, "Elements computed in {} sec", m_data.m_timer_compute);
        Log::info(3, "Elements assembled in {} sec", m_data.m_timer_assemble);
    }

public:     // methods
    Vector hl_inv_v(Ref<const Vector> v)
    {
        if (nb_variables() == 0) {
            return Vector(0);
        }

        if (m_linear_solver.factorize(hl())) {
            throw std::runtime_error("Factorization failed");
        }

        Vector x(nb_variables());

        if (m_linear_solver.solve(v, x)) {
            throw std::runtime_error("Solve failed");
        }

        return x;
    }

    Vector hl_v(Ref<const Vector> v) const
    {
        return hl().selfadjointView<Eigen::Lower>() * v;
    }

    // Pointer<Problem> clone() const
    // {
    //     return new_<Problem>(*this);
    // }

public:     // methods: model properties
    bool parallel() const noexcept
    {
        return m_parallel;
    }

    void set_parallel(const bool value) noexcept
    {
        m_parallel = value;
    }

    bool is_constrained() const noexcept
    {
        return !m_equations.empty();
    }

    index nb_elements_f() const noexcept
    {
        return length(m_elements_f);
    }

    index nb_elements_g() const noexcept
    {
        return length(m_elements_g);
    }

    const std::vector<Pointer<Equation>>& equations() const noexcept
    {
        return m_equations;
    }

    const std::vector<Pointer<Variable>>& variables() const noexcept
    {
        return m_variables;
    }

    index nb_equations() const noexcept
    {
        return length(m_equations);
    }

    index nb_variables() const noexcept
    {
        return length(m_variables);
    }

    const Pointer<Variable>& variable(const index index) const
    {
        return m_variables.at(index);
    }

    const Pointer<Equation>& equation(const index index) const
    {
        return m_equations.at(index);
    }

public:     // methods: input
    Vector delta() const
    {
        Vector result(nb_variables());

        for (index i = 0; i < length(result); i++) {
            result(i) = variable(i)->delta();
        }

        return result;
    }

    void set_delta(Ref<const Vector> value) const
    {
        if (length(value) != nb_variables()) {
            throw std::runtime_error("Invalid size");
        }

        for (index i = 0; i < length(value); i++) {
            variable(i)->set_delta(value[i]);
        }
    }

    void set_delta(double* const value) const
    {
        set_delta(Map<const Vector>(value, nb_variables()));
    }

    Vector x() const
    {
        Vector result(nb_variables());

        for (index i = 0; i < length(result); i++) {
            result(i) = variable(i)->act_value();
        }

        return result;
    }

    void set_x(Ref<const Vector> value) const
    {
        if (length(value) != nb_variables()) {
            throw std::runtime_error("Invalid size");
        }

        for (index i = 0; i < length(value); i++) {
            variable(i)->set_act_value(value[i]);
        }
    }

    void set_x(double* const value) const
    {
        set_x(Map<const Vector>(value, nb_variables()));
    }

    Vector variable_multipliers() const
    {
        Vector result(nb_variables());

        for (index i = 0; i < length(result); i++) {
            result(i) = variable(i)->multiplier();
        }

        return result;
    }

    void set_variable_multipliers(Ref<const Vector> value) const
    {
        if (length(value) != nb_variables()) {
            throw std::runtime_error("Invalid size");
        }

        for (index i = 0; i < length(value); i++) {
            variable(i)->set_multiplier(value[i]);
        }
    }

    void set_variable_multipliers(double* const value) const
    {
        set_variable_multipliers(Map<const Vector>(value, nb_variables()));
    }

    Vector equation_multipliers() const
    {
        Vector result(nb_equations());

        for (index i = 0; i < length(result); i++) {
            result(i) = equation(i)->multiplier();
        }

        return result;
    }

    void set_equation_multipliers(Ref<const Vector> value) const
    {
        if (length(value) != nb_equations()) {
            throw std::runtime_error("Invalid size");
        }

        for (index i = 0; i < length(value); i++) {
            equation(i)->set_multiplier(value[i]);
        }
    }

    void set_equation_multipliers(double* const value) const
    {
        set_equation_multipliers(Map<const Vector>(value, nb_equations()));
    }

    double sigma() const noexcept
    {
        return m_sigma;
    }

    void set_sigma(const double value) noexcept
    {
        m_sigma = value;
    }

public:     // methods: output f
    double f() const noexcept
    {
        return m_data.f();
    }

    void set_f(const double value) noexcept
    {
        m_data.f() = value;
    }

public:     // methods: output g
    Ref<Vector> g() noexcept
    {
        return m_data.g();
    }

    Ref<const Vector> g() const noexcept
    {
        return m_data.g();
    }

    double& g(const index index)
    {
        return m_data.g(index);
    }

    double g(const index index) const
    {
        return m_data.g(index);
    }

public:     // methods: output df
    Ref<Vector> df() noexcept
    {
        return m_data.df();
    }

    Ref<const Vector> df() const noexcept
    {
        return m_data.df();
    }

    double& df(const index index)
    {
        return m_data.df(index);
    }

    double df(const index index) const
    {
        return m_data.df(index);
    }

public:     // methods: output dg
    Map<const Sparse> dg() const noexcept
    {
        return Map<const Sparse>(nb_equations(), nb_variables(), m_dg_structure.nb_nonzeros(), m_dg_structure.ia().data(), m_dg_structure.ja().data(), m_data.dg().data());
    }

    Ref<Vector> dg_values() noexcept
    {
        return m_data.dg();
    }

    Ref<const Vector> dg_values() const noexcept
    {
        return m_data.dg();
    }

    const std::vector<int>& dg_indptr() const noexcept
    {
        return m_dg_structure.ia();
    }

    const std::vector<int>& dg_indices() const noexcept
    {
        return m_dg_structure.ja();
    }

    double& dg(const index index)
    {
        return m_data.dg(index);
    }

    double dg(const index index) const
    {
        return m_data.dg(index);
    }

    double& dg(const index row, const index col)
    {
        const auto index = m_dg_structure.coeff(m_data.dg(), row, col);
        return m_data.dg(index);
    }

    double dg(const index row, const index col) const
    {
        const auto index = m_dg_structure.coeff(m_data.dg(), row, col);
        return m_data.dg(index);
    }

public:     // methods: output hl
    Map<const Sparse> hl() const noexcept
    {
        return Map<const Sparse>(m_hl_structure.rows(), m_hl_structure.cols(), m_hl_structure.nb_nonzeros(), m_hl_structure.ia().data(), m_hl_structure.ja().data(), m_data.hl().data());
    }

    Ref<Vector> hl_values() noexcept
    {
        return m_data.hl();
    }

    Ref<const Vector> hl_values() const noexcept
    {
        return m_data.hl();
    }

    const std::vector<int>& hl_indptr() const noexcept
    {
        return m_hl_structure.ia();
    }

    const std::vector<int>& hl_indices() const noexcept
    {
        return m_hl_structure.ja();
    }

    double& hl(const index index)
    {
        return m_data.hl(index);
    }

    double hl(const index index) const
    {
        return m_data.hl(index);
    }

    double& hl(const index row, const index col)
    {
        const auto index = m_hl_structure.coeff(m_data.hl(), row, col);
        return m_data.hl(index);
    }

    double hl(const index row, const index col) const
    {
        const auto index = m_hl_structure.coeff(m_data.hl(), row, col);
        return m_data.hl(index);
    }

public:     // methods: python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Type = Problem;
        using Holder = Pointer<Type>;

        const std::string name = "Problem";

        py::object scipy_sparse = py::module::import("scipy.sparse");
        py::object csc_matrix = scipy_sparse.attr("csc_matrix");

        py::class_<Type, Holder>(m, name.c_str())
            // constructors
            .def(py::init<ElementsF, ElementsG, Settings>(),
                "objective"_a, "constraints"_a=py::list(),
                "linear_solver"_a=Settings())
            // read-only properties
            .def_property_readonly("is_constrained", &Type::is_constrained)
            .def_property_readonly("equations", &Type::equations)
            .def_property_readonly("variables", &Type::variables)
            .def_property_readonly("g", py::overload_cast<>(&Type::g))
            .def_property_readonly("df", py::overload_cast<>(&Type::df))
            .def_property_readonly("dg", [=](Type& self) {
                return csc_matrix(
                    std::make_tuple(self.dg_values(), self.dg_indices(), self.dg_indptr()),
                    std::make_pair(self.nb_equations(), self.nb_variables())
                ).release();
            })
            .def_property_readonly("dg_values", py::overload_cast<>(&Type::dg_values))
            .def_property_readonly("dg_indptr", &Type::dg_indptr)
            .def_property_readonly("dg_indices", &Type::dg_indices)
            .def_property_readonly("hl", [=](Type& self) {
                return csc_matrix(
                    std::make_tuple(self.hl_values(), self.hl_indices(), self.hl_indptr()),
                    std::make_pair(self.nb_variables(), self.nb_variables())
                ).release();
            })
            .def_property_readonly("hl_values", py::overload_cast<>(&Type::hl_values))
            .def_property_readonly("hl_indptr", &Type::hl_indptr)
            .def_property_readonly("hl_indices", &Type::hl_indices)
            .def_property_readonly("nb_equations", &Type::nb_equations)
            .def_property_readonly("nb_variables", &Type::nb_variables)
            // properties
            .def_property("f", &Type::f, &Type::set_f)
            .def_property("parallel", &Type::parallel, &Type::set_parallel)
            .def_property("sigma", &Type::sigma, &Type::set_sigma)
            .def_property("delta",py::overload_cast<>(&Type::delta, py::const_),
                py::overload_cast<Ref<const Vector>>(&Type::set_delta, py::const_))
            .def_property("x", py::overload_cast<>(&Type::x, py::const_),
                py::overload_cast<Ref<const Vector>>(&Type::set_x, py::const_))
            .def_property("variable_multipliers", py::overload_cast<>(&Type::variable_multipliers, py::const_),
                py::overload_cast<Ref<const Vector>>(&Type::set_variable_multipliers, py::const_))
            .def_property("equation_multipliers", py::overload_cast<>(&Type::equation_multipliers, py::const_),
                py::overload_cast<Ref<const Vector>>(&Type::set_equation_multipliers, py::const_))
            // methods
            // // .def("clone", &Type::clone)
            .def("compute", &Type::compute, "order"_a=2, py::call_guard<py::gil_scoped_release>())
            .def("hl_inv_v", &Type::hl_inv_v)
            .def("hl_v", &Type::hl_v)
        ;
    }
};

} // namespace EQlib