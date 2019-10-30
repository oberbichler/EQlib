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

    int m_nb_threats;
    int m_grainsize;

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

    LinearSolver m_linear_solver;

public:     // constructors
    Problem(ElementsF elements_f, ElementsG elements_g)
    : m_elements_f(std::move(elements_f))
    , m_elements_g(std::move(elements_g))
    , m_sigma(1.0)
    , m_nb_threats(1)
    , m_grainsize(100)
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

        Log::info(2, "The hessian has {} nonzero entries ({:.3f}%)",
            m_hl_structure.nb_nonzeros(), m_hl_structure.density() * 100.0);

        Log::info(2, "The jacobian of the constraints has {} nonzero entries ({:.3f}%)",
            m_dg_structure.nb_nonzeros(), m_dg_structure.density() * 100.0);

        m_data.resize(n, m, m_dg_structure.nb_nonzeros(), m_hl_structure.nb_nonzeros(), m_max_element_n, m_max_element_m);

        Log::info(1, "The problem occupies {} MB", m_data.values().size() * 8.0 / 1'024 / 1'024);


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

            index size_g = TOrder > 0 ? n : 0;
            index size_h = TOrder > 1 ? n : 0;

            Map<Vector> g(data.m_buffer.data(), size_g);
            Map<Matrix> h(data.m_buffer.data() + size_g, size_h, size_h);

            Timer timer_element_compute;

            const double f = element_f.compute(g, h);

            data.computation_time() += timer_element_compute.ellapsed();

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

                    index index = m_hl_structure.get_index(row.global, col.global);

                    data.hl(index) += h(row.local, col.local);
                }
            }

            data.assemble_time() += timer_element_assemble.ellapsed();
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

            Timer timer_element_compute;

            element_g.compute(fs, gs, hs);

            data.computation_time() += timer_element_compute.ellapsed();

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

                    const index dg_value_i = m_dg_structure.get_index(equation_index.global, col.global);

                    data.dg(dg_value_i) += local_g(col.local);

                    if constexpr(TOrder < 2) {
                        continue;
                    }

                    for (index row_i = col_i; row_i < length(variable_indices); row_i++) {
                        const auto row = variable_indices[row_i];

                        const index hl_value_i = m_hl_structure.get_index(row.global, col.global);

                        data.hl(hl_value_i) += local_h(row.local, col.local);
                    }
                }
            }

            data.assemble_time() += timer_element_assemble.ellapsed();
        }
    }

    template <typename TBegin, typename TEnd>
    void compute_elements_f(const index order, ProblemData& data, const TBegin& begin, const TEnd& end)
    {
        switch (order) {
        case 0:
            compute_elements_f<0>(data, begin, end);
            break;
        case 1:
            compute_elements_f<1>(data, begin, end);
            break;
        case 2:
            compute_elements_f<2>(data, begin, end);
            break;
        }
    }

    template <typename TBegin, typename TEnd>
    void compute_elements_g(const index order, ProblemData& data, const TBegin& begin, const TEnd& end)
    {
        switch (order) {
        case 0:
            compute_elements_g<0>(data, begin, end);
            break;
        case 1:
            compute_elements_g<1>(data, begin, end);
            break;
        case 2:
            compute_elements_g<2>(data, begin, end);
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

        tbb::combinable<ProblemData> m_local_data(m_data);

        Log::info(2, "Compute objective...");

        if (m_nb_threats == 1) {
            compute_elements_f(order, m_data, 0, nb_elements_f());
        } else {
            tbb::task_arena arena(m_nb_threats < 1 ? tbb::task_arena::automatic : m_nb_threats);

            arena.execute([&]() {
                tbb::parallel_for(tbb::blocked_range<index>(0, nb_elements_f(), m_grainsize),
                [&](const tbb::blocked_range<index>& range) {
                        auto& data = m_local_data.local();
                        compute_elements_f(order, data, range.begin(), range.end());
                });
            });
        }

        m_data.f() *= sigma();

        if (order > 0) {
            m_data.df() *= sigma();
        }

        if (order > 1) {
            m_data.hl() *= sigma();
        }

        Log::info(2, "Compute constraints...");

        if (m_nb_threats == 1) {
            compute_elements_g(order, m_data, 0, nb_elements_g());
        } else {
            tbb::task_arena arena(m_nb_threats < 1 ? tbb::task_arena::automatic : m_nb_threats);

            arena.execute([&]() {
            tbb::parallel_for(tbb::blocked_range<index>(0, nb_elements_g(), m_grainsize),
                [&](const tbb::blocked_range<index>& range) {
                        auto& data = m_local_data.local();
                        compute_elements_g(order, data, range.begin(), range.end());
                });
            });
        }

        if (m_nb_threats != 1) {
            Log::info(5, "Combine results...");

            m_local_data.combine_each([&](const ProblemData& local) {
                m_data += local;
            });
        }

        Log::info(2, "Problem computed in {} sec", timer.ellapsed());

        Log::info(3, "Element computation took {} sec", m_data.computation_time());
        Log::info(3, "Assembly of the system took {} sec", m_data.assemble_time());
    }

public:     // methods
    Vector hl_inv_v(Ref<const Vector> v)
    {
        if (nb_variables() == 0) {
            return Vector(0);
        }

        Map<const Sparse> hl = this->hl();

        if (m_linear_solver.factorize(hl)) {
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

    void hl_add_diagonal(const double value)
    {
        for (index i = 0; i < nb_variables(); i++) {
            hl(i, i) += value;
        }
    }

    // Pointer<Problem> clone() const
    // {
    //     return new_<Problem>(*this);
    // }

public:     // methods: model properties
    int nb_threats() const noexcept
    {
        return m_nb_threats;
    }

    void set_nb_threats(const int value) noexcept
    {
        m_nb_threats = value;
    }

    int grainsize() const noexcept
    {
        return m_grainsize;
    }

    void set_grainsize(const int value) noexcept
    {
        m_grainsize = value;
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
    Vector x() const
    {
        Vector result(nb_variables());

        for (index i = 0; i < length(result); i++) {
            result(i) = variable(i)->value();
        }

        return result;
    }

    void set_x(Ref<const Vector> value) const
    {
        if (length(value) != nb_variables()) {
            throw std::runtime_error("Invalid size");
        }

        for (index i = 0; i < length(value); i++) {
            variable(i)->set_value(value[i]);
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

public:     // methods: output values
    Ref<Vector> values() noexcept
    {
        return Map<Vector>(m_data.values_ptr(), m_data.values().size());
    }

    Ref<const Vector> values() const noexcept
    {
        return Map<const Vector>(m_data.values_ptr(), m_data.values().size());
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
        const index index = m_dg_structure.get_index(row, col);
        return m_data.dg(index);
    }

    double dg(const index row, const index col) const
    {
        const index index = m_dg_structure.get_index(row, col);
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
        index index = m_hl_structure.get_index(row, col);
        return m_data.hl(index);
    }

    double hl(const index row, const index col) const
    {
        index index = m_hl_structure.get_index(row, col);
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
            .def(py::init<ElementsF, ElementsG>(), "objective"_a=py::list(), "constraints"_a=py::list())
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
            .def_property_readonly("values", py::overload_cast<>(&Type::values))
            // properties
            .def_property("f", &Type::f, &Type::set_f)
            .def_property("nb_threats", &Type::nb_threats, &Type::set_nb_threats)
            .def_property("grainsize", &Type::grainsize, &Type::set_grainsize)
            .def_property("sigma", &Type::sigma, &Type::set_sigma)
            .def_property("x", py::overload_cast<>(&Type::x, py::const_),
                py::overload_cast<Ref<const Vector>>(&Type::set_x, py::const_))
            .def_property("variable_multipliers", py::overload_cast<>(&Type::variable_multipliers, py::const_),
                py::overload_cast<Ref<const Vector>>(&Type::set_variable_multipliers, py::const_))
            .def_property("equation_multipliers", py::overload_cast<>(&Type::equation_multipliers, py::const_),
                py::overload_cast<Ref<const Vector>>(&Type::set_equation_multipliers, py::const_))
            // methods
            // // .def("clone", &Type::clone)
            .def("compute", &Type::compute, "order"_a=2, py::call_guard<py::gil_scoped_release>())
            .def("hl_add_diagonal", &Type::hl_add_diagonal, "values"_a)
            .def("hl_inv_v", &Type::hl_inv_v)
            .def("hl_v", &Type::hl_v)
        ;
    }
};

} // namespace EQlib
