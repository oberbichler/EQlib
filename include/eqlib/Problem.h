#pragma once

#include "Constraint.h"
#include "Define.h"
#include "LinearSolver.h"
#include "Objective.h"
#ifdef EQLIB_USE_MKL
#include "PardisoLDLT.h"
#endif
#include "ProblemData.h"
#include "Settings.h"
#include "SimplicialLDLT.h"
#include "SparseStructure.h"
#include "Timer.h"

#include <sparsehash/dense_hash_map>

#include <tsl/robin_set.h>

#include <set>
#include <tuple>
#include <utility>
#include <vector>

namespace eqlib {

class Problem
{
private:    // types
    using Type = Problem;

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

    int m_nb_threads;
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

    SparseStructure<double, int, true> m_structure_dg;
    SparseStructure<double, int, true> m_structure_hm;

    ProblemData m_data;

    std::unique_ptr<LinearSolver> m_linear_solver;

public:     // constructors
    Problem(ElementsF elements_f, ElementsG elements_g, const index nb_threads=1)
    : m_elements_f(std::move(elements_f))
    , m_elements_g(std::move(elements_g))
    , m_sigma(1.0)
    , m_nb_threads(nb_threads)
    , m_grainsize(100)
    , m_max_element_n(0)
    , m_max_element_m(0)
    {
        Log::task_begin("Initialize problem...");

        Timer timer;

        const auto nb_elements_f = length(m_elements_f);
        const auto nb_elements_g = length(m_elements_g);

        Log::task_info("The objective consists of {} elements", nb_elements_f);
        Log::task_info("The constraints consist of {} elements", nb_elements_g);


        Log::task_step("Getting equations and variables...");

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


        Log::task_step("Creating the set of unique equations...");

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


        Log::task_step("Creating the set of unique variables...");

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

        Log::task_info("The problem contains {} variables", nb_variables);
        Log::task_info("The problem contains {} constraint equations", nb_equations);


        Log::task_step("Compute indices for variables and equations...");

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


        Log::task_step("Compute indices for elements...");

        // variable indices f

        m_element_f_variable_indices.resize(nb_elements_f);
        m_element_g_equation_indices.resize(nb_elements_g);
        m_element_g_variable_indices.resize(nb_elements_g);

        #pragma omp parallel if(m_nb_threads != 1) num_threads(m_nb_threads)
        {
            #pragma omp for schedule(guided, m_grainsize) nowait
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

            #pragma omp for schedule(guided, m_grainsize) nowait
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

            #pragma omp for schedule(guided, m_grainsize)
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
        }


        Log::task_step("Analyse sparse patterns...");

        const auto n = length(m_variables);
        const auto m = length(m_equations);

        std::vector<std::set<index>> m_pattern_dg(m);
        std::vector<std::set<index>> m_pattern_hm(n);

        #pragma omp parallel if(m_nb_threads != 1) num_threads(m_nb_threads)
        {
            std::vector<tsl::robin_set<index>> pattern_dg(m);
            std::vector<tsl::robin_set<index>> pattern_hm(n);

            #pragma omp for schedule(guided, m_grainsize) nowait
            for (index i = 0; i < length(m_elements_f); i++) {
                const auto& variable_indices = m_element_f_variable_indices[i];

                for (index row_i = 0; row_i < length(variable_indices); row_i++) {
                    const auto row = variable_indices[row_i];

                    for (index col_i = row_i; col_i < length(variable_indices); col_i++) {
                        const auto col = variable_indices[col_i];

                        pattern_hm[row.global].insert(col.global);
                    }
                }
            }

            #pragma omp for schedule(guided, m_grainsize) nowait
            for (index i = 0; i < length(m_elements_g); i++) {
                const auto& equation_indices = m_element_g_equation_indices[i];
                const auto& variable_indices = m_element_g_variable_indices[i];

                for (const auto row : equation_indices) {
                    for (const auto col : variable_indices) {
                        pattern_dg[row.global].insert(col.global);
                    }
                }

                for (index row_i = 0; row_i < length(variable_indices); row_i++) {
                    const auto row = variable_indices[row_i];

                    for (index col_i = row_i; col_i < length(variable_indices); col_i++) {
                        const auto col = variable_indices[col_i];

                        pattern_hm[row.global].insert(col.global);
                    }
                }
            }

            #pragma omp critical
            for (index i = 0; i < pattern_hm.size(); i++) {
                for (const auto j : pattern_hm[i]) {
                    m_pattern_hm[i].insert(j);
                }
            }

            #pragma omp critical
            for (index i = 0; i < pattern_dg.size(); i++) {
                for (const auto j : pattern_dg[i]) {
                    m_pattern_dg[i].insert(j);
                }
            }
        }


        Log::task_step("Allocate memory...");

        m_structure_dg.set(m, n, m_pattern_dg);
        m_structure_hm.set(n, n, m_pattern_hm);

        Log::task_info("The hessian has {} nonzero entries ({:.3f}%)",
            m_structure_hm.nb_nonzeros(), m_structure_hm.density() * 100.0);

        Log::task_info("The jacobian of the constraints has {} nonzero entries ({:.3f}%)",
            m_structure_dg.nb_nonzeros(), m_structure_dg.density() * 100.0);

        m_data.resize(n, m, m_structure_dg.nb_nonzeros(), m_structure_hm.nb_nonzeros(), m_max_element_n, m_max_element_m);

        Log::task_info("The problem occupies {} MB", m_data.values().size() * 8.0 / 1'024 / 1'024);

        #ifdef EQLIB_USE_MKL
        m_linear_solver = new_<PardisoLDLT>();
        #else
        m_linear_solver = new_<SimplicialLDLT>();
        #endif

        Log::task_end("Problem initialized in {:.3f} sec", timer.ellapsed());
    }

private:    // methods: computation
    template <index TOrder>
    void compute_element_f(ProblemData& data, const index i)
    {
        static_assert(0 <= TOrder && TOrder <= 2);

        const auto& element_f = *m_elements_f[i];

        if (!element_f.is_active()) {
            return;
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
            return;
        }

            for (index row_i = 0; row_i < length(variable_indices); row_i++) {
                const auto row = variable_indices[row_i];

                data.df(row.global) += g(row.local);

            if constexpr(TOrder < 2) {
                return;
            }

                for (index col_i = row_i; col_i != length(variable_indices); ++col_i) {
                    const auto col = variable_indices[col_i];

                index index = m_structure_hm.get_index(row.global, col.global);

                data.hm_value(index) += h(row.local, col.local);
            }
        }

        data.assemble_time() += timer_element_assemble.ellapsed();
    }

    template <index TOrder>
    void compute_element_g(ProblemData& data, const index i)
    {
        static_assert(0 <= TOrder && TOrder <= 2);

        const auto& element_g = *m_elements_g[i];

        if (!element_g.is_active()) {
            return;
        }

        const auto& equation_indices = m_element_g_equation_indices[i];
        const auto& variable_indices = m_element_g_variable_indices[i];

        const auto m = m_element_g_nb_equations[i];
        const auto n = m_element_g_nb_variables[i];

        if (n == 0 || m == 0) { // FIXME: check reduces counts
            return;
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
                return;
            }

            auto& local_g = gs[equation_index.local];
            auto& local_h = hs[equation_index.local];

            local_h *= equation->multiplier();

                for (index row_i = 0; row_i < length(variable_indices); row_i++) {
                    const auto row = variable_indices[row_i];

                    const index dg_value_i = m_structure_dg.get_index(equation_index.global, row.global);

                    data.dg_value(dg_value_i) += local_g(row.local);

                if constexpr(TOrder < 2) {
                    return;
                }

                    for (index col_i = row_i; col_i < length(variable_indices); col_i++) {
                        const auto col = variable_indices[col_i];

                    const index hm_value_i = m_structure_hm.get_index(row.global, col.global);

                    data.hm_value(hm_value_i) += local_h(row.local, col.local);
                }
            }
        }

        data.assemble_time() += timer_element_assemble.ellapsed();
    }

public:     // methods: computation
    template <bool TParallel, bool TInfo, index TOrder>
    void compute()
    {
        static_assert(0 <= TOrder && TOrder <= 2);

        if constexpr(TInfo) {
            Log::task_begin("Compute problem...");
        }

        Timer timer;

        m_data.set_zero();

        if constexpr(TParallel) {
            ProblemData local_data(m_data);

            // FIXME: use private(m_data) when msvc supports omp private with member variables
            #pragma omp parallel if(m_nb_threads != 1) num_threads(m_nb_threads) firstprivate(local_data)
            {
                #pragma omp for schedule(guided, m_grainsize) nowait
                for (index i = 0; i < nb_elements_f(); i++) {
                    compute_element_f<TOrder>(local_data, i);
                }

                if (sigma() != 1.0) {
                    local_data.f() *= sigma();

                    if constexpr(TOrder > 0) {
                        local_data.df() *= sigma();
                    }

                    if constexpr(TOrder > 1) {
                        local_data.hm_values() *= sigma();
                    }
                }

                #pragma omp for schedule(guided, m_grainsize) nowait
                for (index i = 0; i < nb_elements_g(); i++) {
                    compute_element_g<TOrder>(local_data, i);
                }

                #pragma omp critical
                m_data += local_data;
            }
        } else {
            for (index i = 0; i < nb_elements_f(); i++) {
                compute_element_f<TOrder>(m_data, i);
            }

            if (sigma() != 1.0) {
                m_data.f() *= sigma();

                if constexpr(TOrder > 0) {
                    m_data.df() *= sigma();
                }

                if constexpr(TOrder > 1) {
                    m_data.hm_values() *= sigma();
                }
            }

            for (index i = 0; i < nb_elements_g(); i++) {
                compute_element_g<TOrder>(m_data, i);
            }
        }

        if constexpr(TInfo) {
            Log::task_info("Element computation took {} sec", m_data.computation_time());
            Log::task_info("Assembly of the system took {} sec", m_data.assemble_time());

            Log::task_end("Problem computed in {:.3f} sec", timer.ellapsed());
        }
    }

    template <bool TInfo, index TOrder>
    void compute()
    {
        if (m_nb_threads == 1) {
            compute<false, TInfo, TOrder>();
        } else {
            compute<true, TInfo, TOrder>();
        }
    }

    template <bool TInfo>
    void compute(const index order = 2)
    {
        if (m_nb_threads == 1) {
            switch (order) {
            case 0:
                compute<false, TInfo, 0>();
                break;
            case 1:
                compute<false, TInfo, 1>();
                break;
            case 2:
                compute<false, TInfo, 2>();
                break;
            default:
                throw std::invalid_argument("order");
            }
        } else {
            switch (order) {
            case 0:
                compute<true, TInfo, 0>();
                break;
            case 1:
                compute<true, TInfo, 1>();
                break;
            case 2:
                compute<true, TInfo, 2>();
                break;
            default:
                throw std::invalid_argument("order");
            }
        }
    }

    void compute(const index order = 2)
    {
        if (m_nb_threads == 1) {
            switch (order) {
            case 0:
                compute<false, true, 0>();
                break;
            case 1:
                compute<false, true, 1>();
                break;
            case 2:
                compute<false, true, 2>();
                break;
            default:
                throw std::invalid_argument("order");
            }
        } else {
            switch (order) {
            case 0:
                compute<true, true, 0>();
                break;
            case 1:
                compute<true, true, 1>();
                break;
            case 2:
                compute<true, true, 2>();
                break;
            default:
                throw std::invalid_argument("order");
            }
        }
    }

public:     // methods
    Vector hm_inv_v(Ref<const Vector> v)
    {
        if (nb_variables() == 0) {
            return Vector(0);
        }

        if (m_linear_solver->factorize(m_structure_hm.ia(), m_structure_hm.ja(), m_data.hm_values())) {
            throw std::runtime_error("Factorization failed");
        }

        Vector x(nb_variables());

        if (m_linear_solver->solve(m_structure_hm.ia(), m_structure_hm.ja(), m_data.hm_values(), v, x)) {
            throw std::runtime_error("Solve failed");
        }

        return x;
    }

    Vector hm_v(Ref<const Vector> v) const
    {
        return hm().selfadjointView<Eigen::Upper>() * v;
    }

    void hm_add_diagonal(const double value)
    {
        for (index i = 0; i < nb_variables(); i++) {
            hm(i, i) += value;
        }
    }

    // Pointer<Problem> clone() const
    // {
    //     return new_<Problem>(*this);
    // }

    std::string solver_name() const
    {
        return m_linear_solver->solver_name();
    }

public:     // methods: model properties
    int nb_threads() const noexcept
    {
        return m_nb_threads;
    }

    void set_nb_threads(const int value) noexcept
    {
        m_nb_threads = value;
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

    std::vector<std::pair<double, double>> equation_bounds() const {
        std::vector<std::pair<double, double>> bounds(nb_equations());

        for (index i = 0; i < nb_equations(); i++) {
            bounds[i] = {equation(i)->lower_bound(), equation(i)->upper_bound()};
        }

        return bounds;
    }

    std::vector<std::pair<double, double>> variable_bounds() const {
        std::vector<std::pair<double, double>> bounds(nb_variables());

        for (index i = 0; i < nb_variables(); i++) {
            bounds[i] = {variable(i)->lower_bound(), variable(i)->upper_bound()};
        }

        return bounds;
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
    Ref<const Sparse> dg() const noexcept
    {
        return Map<const Sparse>(nb_equations(), nb_variables(), m_structure_dg.nb_nonzeros(), m_structure_dg.ia().data(), m_structure_dg.ja().data(), m_data.dg_ptr());
    }

    Ref<Vector> dg_values() noexcept
    {
        return m_data.dg_values();
    }

    Ref<const Vector> dg_values() const noexcept
    {
        return m_data.dg_values();
    }

    const std::vector<int>& dg_indptr() const noexcept
    {
        return m_structure_dg.ia();
    }

    const std::vector<int>& dg_indices() const noexcept
    {
        return m_structure_dg.ja();
    }

    double& dg(const index index)
    {
        return m_data.dg_value(index);
    }

    double dg(const index index) const
    {
        return m_data.dg_value(index);
    }

    double& dg(const index row, const index col)
    {
        const index index = m_structure_dg.get_index(row, col);
        return m_data.dg_value(index);
    }

    double dg(const index row, const index col) const
    {
        const index index = m_structure_dg.get_index(row, col);
        return m_data.dg_value(index);
    }

public:     // methods: output hm
    Map<const Sparse> hm() const noexcept
    {
        return Map<const Sparse>(m_structure_hm.rows(), m_structure_hm.cols(), m_structure_hm.nb_nonzeros(), m_structure_hm.ia().data(), m_structure_hm.ja().data(), m_data.hm_ptr());
    }

    Ref<Vector> hm_values() noexcept
    {
        return m_data.hm_values();
    }

    Ref<const Vector> hm_values() const noexcept
    {
        return m_data.hm_values();
    }

    const std::vector<int>& hm_indptr() const noexcept
    {
        return m_structure_hm.ia();
    }

    const std::vector<int>& hm_indices() const noexcept
    {
        return m_structure_hm.ja();
    }

    double& hm(const index index)
    {
        return m_data.hm_value(index);
    }

    double hm(const index index) const
    {
        return m_data.hm_value(index);
    }

    double& hm(const index row, const index col)
    {
        index index = m_structure_hm.get_index(row, col);
        return m_data.hm_value(index);
    }

    double hm(const index row, const index col) const
    {
        index index = m_structure_hm.get_index(row, col);
        return m_data.hm_value(index);
    }

public:     // methods: python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Holder = Pointer<Type>;

        const std::string name = "Problem";

        py::object scipy_sparse = py::module::import("scipy.sparse");
        py::object csr_matrix = scipy_sparse.attr("csr_matrix");

        py::class_<Type, Holder>(m, name.c_str())
            // constructors
            .def(py::init<ElementsF, ElementsG, index>(), "objective"_a=py::list(), "constraints"_a=py::list(),
                "nb_threads"_a=1)
            // read-only properties
            .def_property_readonly("is_constrained", &Type::is_constrained)
            .def_property_readonly("equations", &Type::equations)
            .def_property_readonly("variables", &Type::variables)
            .def_property_readonly("g", py::overload_cast<>(&Type::g))
            .def_property_readonly("df", py::overload_cast<>(&Type::df))
            .def_property_readonly("dg", [=](Type& self) {
                return csr_matrix(
                    std::make_tuple(self.dg_values(), self.dg_indices(), self.dg_indptr()),
                    std::make_pair(self.nb_equations(), self.nb_variables())
                ).release();
            })
            .def_property_readonly("dg_values", py::overload_cast<>(&Type::dg_values))
            .def_property_readonly("dg_indptr", &Type::dg_indptr)
            .def_property_readonly("dg_indices", &Type::dg_indices)
            .def_property_readonly("hm", [=](Type& self) {
                return csr_matrix(
                    std::make_tuple(self.hm_values(), self.hm_indices(), self.hm_indptr()),
                    std::make_pair(self.nb_variables(), self.nb_variables())
                ).release();
            })
            .def_property_readonly("hm_values", py::overload_cast<>(&Type::hm_values))
            .def_property_readonly("hm_indptr", &Type::hm_indptr)
            .def_property_readonly("hm_indices", &Type::hm_indices)
            .def_property_readonly("nb_equations", &Type::nb_equations)
            .def_property_readonly("nb_variables", &Type::nb_variables)
            .def_property_readonly("values", py::overload_cast<>(&Type::values))
            .def_property_readonly("equation_bounds", &Type::equation_bounds)
            .def_property_readonly("variable_bounds", &Type::variable_bounds)
            // properties
            .def_property("f", &Type::f, &Type::set_f)
            .def_property("nb_threads", &Type::nb_threads, &Type::set_nb_threads)
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
            .def("compute", &Type::compute<true>, "order"_a=2, py::call_guard<py::gil_scoped_release>())
            .def("hm_add_diagonal", &Type::hm_add_diagonal, "value"_a)
            .def("hm_inv_v", &Type::hm_inv_v)
            .def("hm_v", &Type::hm_v)
            .def("f_of", [](Type& self, Ref<const Vector> x) -> double {
                self.set_x(x);
                self.compute<false>(0);
                return self.f();
            }, "x"_a)
            .def("g_of", [](Type& self, Ref<const Vector> x) -> Vector {
                self.set_x(x);
                self.compute<false>(0);
                return self.g();
            }, "x"_a)
            .def("df_of", [](Type& self, Ref<const Vector> x) -> Vector {
                self.set_x(x);
                self.compute<false>(1);
                return self.df();
            }, "x"_a)
            .def("dg_of", [=](Type& self, Ref<const Vector> x) {
                self.set_x(x);
                self.compute<false>(1);
                return csr_matrix(
                    std::make_tuple(self.dg_values(), self.dg_indices(), self.dg_indptr()),
                    std::make_pair(self.nb_equations(), self.nb_variables())
                ).release();
            }, "x"_a)
            .def("hm_of", [=](Type& self, Ref<const Vector> x) {
                self.set_x(x);
                self.compute<false>(2);
                return csr_matrix(
                    std::make_tuple(self.hm_values(), self.hm_indices(), self.hm_indptr()),
                    std::make_pair(self.nb_variables(), self.nb_variables())
                ).release();
            }, "x"_a)
        ;
    }
};

} // namespace eqlib
