#pragma once

#include "common.h"
#include "constraint.h"
#include "objective.h"
#include "request.h"
#include "sparse_structure.h"
#include "parameter.h"
#ifdef EQLIB_USE_MKL
#include "linear_solvers/pardiso_ldlt.h"
#else
#include "linear_solvers/simplicial_ldlt.h"
#endif

#include <algorithm>
#include <future>
#include <mutex>
#include <stdexcept>

namespace eqlib {

struct Index {
    index local;
    index global;

    bool operator<(const Index& other) const noexcept
    {
        return global < other.global;
    }
};

struct ObjectiveIndices {
    ObjectiveIndices() = default;

    ObjectiveIndices(const Pointer<Objective>& objective, const index nb_variables, const index nb_constants)
        : m_objective(objective)
        , variable_indices(nb_variables)
        , constant_indices(nb_constants)
    {
    }

    // objective

    Pointer<Objective> m_objective;

    Objective& objective()
    {
        return *m_objective;
    }

    // variable indices

    std::vector<Index> variable_indices;

    // constant indices

    std::vector<Index> constant_indices;

    // sort

    void sort()
    {
        std::sort(variable_indices.begin(), variable_indices.end());
        std::sort(constant_indices.begin(), constant_indices.end());
    }
};

struct ConstraintIndices {
    ConstraintIndices() = default;

    ConstraintIndices(const Pointer<Constraint>& constraint, const index nb_equations, const index nb_variables, const index nb_constants)
        : m_constraint(constraint)
        , equation_indices(nb_equations)
        , variable_indices(nb_variables)
        , constant_indices(nb_constants)
    {
    }

    // constraint

    Pointer<Constraint> m_constraint;

    Constraint& constraint()
    {
        return *m_constraint;
    }

    // equation indices

    std::vector<Index> equation_indices;

    // variable indices

    std::vector<Index> variable_indices;

    // constant indices

    std::vector<Index> constant_indices;

    // sort

    void sort()
    {
        std::sort(equation_indices.begin(), equation_indices.end());
        std::sort(variable_indices.begin(), variable_indices.end());
        std::sort(constant_indices.begin(), constant_indices.end());
    }
};

struct Problem {
    Problem()
        : m_buffer_size(0)
        , m_nb_nonzeros_dg(0)
        , m_nb_nonzeros_hm(0)
        , m_parallel_tasks(0)
        , m_parallel_block_size(256)
        , m_sigma(1.0)
    {
    }

    static Problem create(const std::vector<Pointer<Objective>>& objectives, const std::vector<Pointer<Constraint>>& constraints)
    {
        // Logger::info("Create a new system");

        Problem result;

        const index nb_objectives = len(objectives);
        const index nb_constraints = len(constraints);

        // Logger::info("{} objective elements", nb_objectives);
        // Logger::info("{} constraints elements", nb_constraints);

        result.m_objective_indices = std::vector<ObjectiveIndices>(nb_objectives);
        result.m_constraint_indices = std::vector<ConstraintIndices>(nb_constraints);

        // index parameters and equations

        for (index i = 0; i < nb_objectives; i++) {
            const auto& objective = objectives[i];

            index local_nb_constants = 0;

            for (const auto& parameter : objective->parameters()) {
                if (parameter->is_active()) {
                    result.m_variable_indices.try_emplace(parameter, len(result.m_variable_indices));
                } else {
                    local_nb_constants += 1;
                    result.m_constant_indices.try_emplace(parameter, len(result.m_constant_indices));
                }
            }

            const index local_nb_variables = objective->nb_parameters() - local_nb_constants;

            result.m_objective_indices[i] = ObjectiveIndices(objective, local_nb_variables, local_nb_constants);
        }

        for (index i = 0; i < nb_constraints; i++) {
            const auto& constraint = constraints[i];

            index local_nb_equations = 0;

            for (const auto& equation : constraint->equations()) {
                if (equation->is_active()) {
                    local_nb_equations += 1;
                    result.m_equation_indices.try_emplace(equation, len(result.m_equation_indices));
                }
            }

            index local_nb_constants = 0;

            for (const auto& parameter : constraint->parameters()) {
                if (parameter->is_active()) {
                    result.m_variable_indices.try_emplace(parameter, len(result.m_variable_indices));
                } else {
                    local_nb_constants += 1;
                    result.m_constant_indices.try_emplace(parameter, len(result.m_constant_indices));
                }
            }

            const index local_nb_variables = constraint->nb_parameters() - local_nb_constants;

            result.m_constraint_indices[i] = ConstraintIndices(constraint, local_nb_equations, local_nb_variables, local_nb_constants);
        }

        // collect equations

        index nb_equations = len(result.m_equation_indices);

        // Logger::info("{} equations", nb_equations);

        result.m_equations = Equations(nb_equations);

        for (const auto& [equation, index] : result.m_equation_indices) {
            result.m_equations[index] = equation;
        }

        // collect variables

        index nb_variables = len(result.m_variable_indices);

        // Logger::info("{} variables", nb_variables);

        result.m_variables = Parameters(nb_variables);

        for (const auto& [variable, index] : result.m_variable_indices) {
            result.m_variables[index] = variable;
        }

        // collect constants

        index nb_constants = len(result.m_constant_indices);

        result.m_constants = Parameters(nb_constants);

        for (const auto& [constant, index] : result.m_constant_indices) {
            result.m_constants[index] = constant;
        }

        // assign indices

        for (auto& objective_indices : result.m_objective_indices) {
            const auto& objective = objective_indices.objective();

            auto variable_indices_it = objective_indices.variable_indices.begin();
            auto constant_indices_it = objective_indices.constant_indices.begin();

            for (index local_index = 0; local_index < objective.nb_parameters(); local_index++) {
                const auto parameter = objective.parameter(local_index);

                if (parameter->is_active()) {
                    const index global_index = result.m_variable_indices[parameter];
                    *(variable_indices_it++) = { local_index, global_index };
                } else {
                    const index global_index = result.m_constant_indices[parameter];
                    *(constant_indices_it++) = { local_index, global_index };
                }
            }

            objective_indices.sort();
        }

        for (auto& constraint_indices : result.m_constraint_indices) {
            const auto& constraint = constraint_indices.constraint();

            auto equation_indices_it = constraint_indices.equation_indices.begin();

            for (index local_index = 0; local_index < constraint.nb_equations(); local_index++) {
                const auto equation = constraint.equation(local_index);

                if (equation->is_active()) {
                    const index global_index = result.m_equation_indices[equation];
                    *(equation_indices_it++) = { local_index, global_index };
                }
            }

            auto variable_indices_it = constraint_indices.variable_indices.begin();
            auto constant_indices_it = constraint_indices.constant_indices.begin();

            for (index local_index = 0; local_index < constraint.nb_parameters(); local_index++) {
                const auto parameter = constraint.parameter(local_index);

                if (parameter->is_active()) {
                    const index global_index = result.m_variable_indices[parameter];
                    *(variable_indices_it++) = { local_index, global_index };
                } else {
                    const index global_index = result.m_constant_indices[parameter];
                    *(constant_indices_it++) = { local_index, global_index };
                }
            }

            constraint_indices.sort();
        }

        // pattern

        std::vector<RobinSet<index>> pattern_dg(nb_equations);
        std::vector<RobinSet<index>> pattern_hm(nb_variables);

        for (auto& constraint_indices : result.m_constraint_indices) {
            const auto& equation_indices = constraint_indices.equation_indices;
            const auto& variable_indices = constraint_indices.variable_indices;

            for (auto& index_i : equation_indices) {
                for (auto& index_j : variable_indices) {
                    pattern_dg[index_i.global].insert(index_j.global);
                }
            }

            const index nb_variables = len(variable_indices);

            for (index i = 0; i < nb_variables; i++) {
                const index global_i = variable_indices[i].global;

                for (index j = i; j < nb_variables; j++) {
                    const index global_j = variable_indices[j].global;

                    pattern_hm[global_i].insert(global_j);
                }
            }

            // update buffer size

            const auto& constraint = constraint_indices.constraint();

            const index required_buffer_size = constraint.nb_equations() + constraint.nb_equations() * constraint.nb_parameters() + constraint.nb_parameters() * constraint.nb_parameters();

            if (required_buffer_size > result.m_buffer_size) {
                result.m_buffer_size = required_buffer_size;
            }
        }

        for (auto& objective_indices : result.m_objective_indices) {
            const index nb_variables = len(objective_indices.variable_indices);

            // update pattern

            for (index i = 0; i < nb_variables; i++) {
                const index global_i = objective_indices.variable_indices[i].global;

                for (index j = i; j < nb_variables; j++) {
                    const index global_j = objective_indices.variable_indices[j].global;

                    pattern_hm[global_i].insert(global_j);
                }
            }

            // update buffer size

            const auto& objective = objective_indices.objective();

            const index required_buffer_size = objective.nb_parameters() + objective.nb_parameters() * objective.nb_parameters();

            if (required_buffer_size > result.m_buffer_size) {
                result.m_buffer_size = required_buffer_size;
            }
        }

        for (const auto& i : pattern_dg) {
            result.m_nb_nonzeros_dg += len(i);
        }

        for (const auto& i : pattern_hm) {
            result.m_nb_nonzeros_hm += len(i);
        }

        result.m_structure_dg = CsrStructure::from_pattern(nb_equations, nb_variables, pattern_dg);
        result.m_structure_hm = CsrStructure::from_pattern(nb_variables, nb_variables, pattern_hm);

        return result;
    }

    static Problem like(const Problem& prototype, const std::vector<Pointer<Objective>>& objectives, const std::vector<Pointer<Constraint>>& constraints)
    {
        Problem result;

        result.m_equation_indices = prototype.m_equation_indices;
        result.m_variable_indices = prototype.m_variable_indices;
        result.m_constant_indices = prototype.m_constant_indices;

        result.m_equations = prototype.m_equations;
        result.m_variables = prototype.m_variables;
        result.m_constants = prototype.m_constants;

        result.m_buffer_size = prototype.m_buffer_size;
        result.m_nb_nonzeros_dg = prototype.m_nb_nonzeros_dg;
        result.m_nb_nonzeros_hm = prototype.m_nb_nonzeros_hm;

        const index nb_objectives = len(objectives);
        const index nb_constraints = len(constraints);

        result.m_objective_indices = std::vector<ObjectiveIndices>(nb_objectives);
        result.m_constraint_indices = std::vector<ConstraintIndices>(nb_constraints);

        // index parameters and equations

        for (index i = 0; i < nb_objectives; i++) {
            const auto& objective = objectives[i];

            index local_nb_constants = 0;

            for (const auto& variable : objective->parameters()) {
                if (!variable->is_active()) {
                    local_nb_constants += 1;
                }
            }

            const index local_nb_variables = objective->nb_parameters() - local_nb_constants;

            result.m_objective_indices[i] = ObjectiveIndices(objective, local_nb_variables, local_nb_constants);
        }

        for (index i = 0; i < nb_constraints; i++) {
            const auto& constraint = constraints[i];

            index local_nb_equations = 0;

            for (const auto& equation : constraint->equations()) {
                if (equation->is_active()) {
                    local_nb_equations += 1;
                }
            }

            index local_nb_constants = 0;

            for (const auto& variable : constraint->parameters()) {
                if (!variable->is_active()) {
                    local_nb_constants += 1;
                }
            }

            const index local_nb_variables = constraint->nb_parameters() - local_nb_constants;

            result.m_constraint_indices[i] = ConstraintIndices(constraint, local_nb_equations, local_nb_variables, local_nb_constants);
        }

        // assign indices

        for (auto& objective_indices : result.m_objective_indices) {
            const auto& objective = objective_indices.objective();

            auto variable_indices_it = objective_indices.variable_indices.begin();
            auto constant_indices_it = objective_indices.constant_indices.begin();

            for (index local_index = 0; local_index < objective.nb_parameters(); local_index++) {
                const auto parameter = objective.parameter(local_index);

                if (parameter->is_active()) {
                    const auto it = prototype.m_variable_indices.find(parameter);

                    if (it != prototype.m_variable_indices.end()) {
                        const index global_index = it.value();
                        *(variable_indices_it++) = { local_index, global_index };
                    }
                } else {
                    const auto it = prototype.m_constant_indices.find(parameter);

                    if (it != prototype.m_constant_indices.end()) {
                        const index global_index = it.value();
                        *(constant_indices_it++) = { local_index, global_index };
                    }
                }
            }

            objective_indices.sort();
        }

        for (auto& constraint_indices : result.m_constraint_indices) {
            const auto& constraint = constraint_indices.constraint();

            auto equation_indices_it = constraint_indices.equation_indices.begin();

            for (index local_index = 0; local_index < constraint.nb_equations(); local_index++) {
                const auto equation = constraint.equation(local_index);

                if (equation->is_active()) {
                    const auto it = prototype.m_equation_indices.find(equation);

                    if (it != prototype.m_equation_indices.end()) {
                        const index global_index = it.value();
                        *(equation_indices_it++) = { local_index, global_index };
                    }
                }
            }

            auto variable_indices_it = constraint_indices.variable_indices.begin();
            auto constant_indices_it = constraint_indices.constant_indices.begin();

            for (index local_index = 0; local_index < constraint.nb_parameters(); local_index++) {
                const auto parameter = constraint.parameter(local_index);

                if (parameter->is_active()) {
                    const auto it = prototype.m_variable_indices.find(parameter);

                    if (it != prototype.m_variable_indices.end()) {
                        const index global_index = it.value();
                        *(variable_indices_it++) = { local_index, global_index };
                    }
                } else {
                    const auto it = prototype.m_constant_indices.find(parameter);

                    if (it != prototype.m_constant_indices.end()) {
                        const index global_index = it.value();
                        *(constant_indices_it++) = { local_index, global_index };
                    }
                }
            }

            constraint_indices.sort();
        }

        return result;
    }

    // constants

    const inline static Map<Vector> EmptyVector = Map<Vector>(nullptr, 0);

    // variable indices

    RobinMap<Pointer<Parameter>, index> m_variable_indices;

    index variable_index(const Pointer<Parameter>& variable) const
    {
        const auto it = m_variable_indices.find(variable);

        if (it == m_variable_indices.end())
            return -1;

        return it.value();
    }

    // constant indices

    RobinMap<Pointer<Parameter>, index> m_constant_indices;

    index constant_index(const Pointer<Parameter>& variable) const
    {
        const auto it = m_constant_indices.find(variable);

        if (it == m_constant_indices.end())
            return -1;

        return it.value();
    }

    // equation indices

    RobinMap<Pointer<Equation>, index> m_equation_indices;

    index equation_index(const Pointer<Equation>& equation) const
    {
        const auto it = m_equation_indices.find(equation);

        if (it == m_equation_indices.end())
            return -1;

        return it.value();
    }

    // objective indices

    std::vector<ObjectiveIndices> m_objective_indices;

    index nb_objectives() const
    {
        return len(m_objective_indices);
    }

    // constraint indices

    std::vector<ConstraintIndices> m_constraint_indices;

    index nb_constraints() const
    {
        return len(m_constraint_indices);
    }

    // equations

    Equations m_equations;

    index nb_equations() const
    {
        return len(m_equations);
    }

    const Equations& equations() const
    {
        return m_equations;
    }

    // variables

    Parameters m_variables;

    index nb_variables() const
    {
        return len(m_variables);
    }

    const Parameters& variables() const
    {
        return m_variables;
    }

    // constants

    Parameters m_constants;

    index nb_constants() const
    {
        return len(m_constants);
    }

    const Parameters& constants() const
    {
        return m_constants;
    }

    // structure_dg

    CsrStructure m_structure_dg;

    CsrStructure structure_dg() const
    {
        return m_structure_dg;
    }

    // structure_hm

    CsrStructure m_structure_hm;

    CsrStructure structure_hm() const
    {
        return m_structure_hm;
    }

    // sigma

    double m_sigma;

    double sigma() const
    {
        return m_sigma;
    }

    // compute f

    template <Request R>
    void _compute_block_f(const index block_begin, const index block_end, Ref<Vector> buffer, double& f, Ref<Vector> df, Ref<Vector> hm)
    {
        constexpr bool request_f = (R & Request::F) != 0;
        constexpr bool request_df = (R & Request::DF) != 0;
        constexpr bool request_hf = (R & Request::HF) != 0;

        for (index i = block_begin; i < block_end; i++) {
            auto& objective_indices = m_objective_indices[i];

            auto& objective = objective_indices.objective();

            if (!objective.is_active())
                continue;

            const index nb_parameters = objective.nb_parameters();

            Map<Vector> local_df(buffer.data(), nb_parameters);
            Map<Matrix> local_hm(buffer.data() + nb_parameters, nb_parameters, nb_parameters);

            const double local_f = objective.compute(local_df, local_hm, R);

            if constexpr (request_f) {
                f += local_f;
            }

            if constexpr (request_df) {
                for (index i = 0; i < len(objective_indices.variable_indices); i++) {
                    const auto index_i = objective_indices.variable_indices[i];

                    df(index_i.global) += local_df(index_i.local);
                }
            }

            if constexpr (request_hf) {
                for (index i = 0; i < len(objective_indices.variable_indices); i++) {
                    const auto index_i = objective_indices.variable_indices[i];

                    for (index j = i; j < len(objective_indices.variable_indices); j++) {
                        const auto index_j = objective_indices.variable_indices[j];

                        const index idx = m_structure_hm.get_index(index_i.global, index_j.global);

                        assert(0 <= idx && idx < len(hm));

                        hm(idx) += local_hm(index_i.local, index_j.local);
                    }
                }
            }
        }
    }

    void _compute_block_f(Request request, const index block_begin, const index block_end, Ref<Vector> buffer, double& f, Ref<Vector> df, Ref<Vector> hm)
    {
        switch (request & 0b00001111) {
        case 0b00000001:
            _compute_block_f<static_cast<Request>(0b00000001)>(block_begin, block_end, buffer, f, df, hm);
            break;
        case 0b00000010:
            _compute_block_f<static_cast<Request>(0b00000010)>(block_begin, block_end, buffer, f, df, hm);
            break;
        case 0b00000011:
            _compute_block_f<static_cast<Request>(0b00000011)>(block_begin, block_end, buffer, f, df, hm);
            break;
        case 0b00000100:
            _compute_block_f<static_cast<Request>(0b00000100)>(block_begin, block_end, buffer, f, df, hm);
            break;
        case 0b00000101:
            _compute_block_f<static_cast<Request>(0b00000101)>(block_begin, block_end, buffer, f, df, hm);
            break;
        case 0b00000110:
            _compute_block_f<static_cast<Request>(0b00000110)>(block_begin, block_end, buffer, f, df, hm);
            break;
        case 0b00000111:
            _compute_block_f<static_cast<Request>(0b00000111)>(block_begin, block_end, buffer, f, df, hm);
            break;
        }
    }

    void _compute_f_async(std::mutex& mutex, Request request, std::atomic<index>& begin, const index end, double& f, Ref<Vector> df, Ref<Vector> hm)
    {
        Vector future_buffer;
        double future_f = 0.0;
        Vector future_df;
        Vector future_hm;

        if (request & Request::DF)
            future_df = Vector::Zero(len(df));

        if (request & Request::HF)
            future_hm = Vector::Zero(len(hm));

        if ((request & Request::DF) || (request & Request::HF))
            future_buffer.resize(buffer_size());

        while (true) {
            index block_begin = begin.fetch_add(m_parallel_block_size);
            index block_end = std::min(block_begin + m_parallel_block_size, end);

            if (block_begin >= end)
                break;

            _compute_block_f(request, block_begin, block_end, future_buffer, future_f, future_df, future_hm);
        }

        mutex.lock();

        if (request & Request::F)
            f += future_f;

        if (request & Request::DF)
            df += future_df;

        if (request & Request::HF)
            hm += future_hm;

        mutex.unlock();
    }

    // compute g

    template <Request R>
    void _compute_block_g(const index block_begin, const index block_end, Ref<Vector> buffer, Ref<Vector> g, Ref<Vector> dg, Ref<Vector> hm)
    {
        constexpr bool request_g = (R & Request::G) != 0;
        constexpr bool request_dg = (R & Request::DG) != 0;
        constexpr bool request_hg = (R & Request::HG) != 0;

        if constexpr (!request_g && !request_dg && !request_hg)
            return;

        for (index i = block_begin; i < block_end; i++) {
            auto& constraint_indices = m_constraint_indices[i];

            auto& constraint = constraint_indices.constraint();

            if (!constraint.is_active())
                continue;

            const index nb_equations = constraint.nb_equations();
            const index nb_parameters = constraint.nb_parameters();

            Map<Vector> local_g(buffer.data(), nb_equations);
            Map<Matrix> local_dg(buffer.data() + nb_equations, nb_equations, nb_parameters);
            Map<Matrix> local_hm(buffer.data() + nb_equations + nb_equations * nb_parameters, nb_parameters, nb_parameters);

            constraint.compute(local_g, local_dg, local_hm, R);

            if constexpr (request_g) {
                for (index i = 0; i < len(constraint_indices.equation_indices); i++) {
                    const auto index_i = constraint_indices.equation_indices[i];

                    if constexpr (request_g)
                        g(index_i.global) += local_g(index_i.local);

                    if constexpr (request_dg) {
                        for (index j = 0; j < len(constraint_indices.variable_indices); j++) {
                            const auto index_j = constraint_indices.variable_indices[j];

                            const index idx = m_structure_dg.get_index(index_i.global, index_j.global);

                            assert(idx >= 0);

                            dg(idx) += local_dg(index_i.local, index_j.local);
                        }
                    }
                }

                if constexpr (request_hg) {
                    for (index i = 0; i < len(constraint_indices.variable_indices); i++) {
                        const auto index_i = constraint_indices.variable_indices[i];

                        for (index j = i; j < len(constraint_indices.variable_indices); j++) {
                            const auto index_j = constraint_indices.variable_indices[j];

                            const index idx = m_structure_hm.get_index(index_i.global, index_j.global);

                            assert(0 <= idx && idx < len(hm));

                            hm(idx) += local_hm(index_i.local, index_j.local);
                        }
                    }
                }
            }
        }
    }

    void _compute_block_g(Request request, const index block_begin, const index block_end, Ref<Vector> buffer, Ref<Vector> g, Ref<Vector> dg, Ref<Vector> hm)
    {
        switch (request & 0b11110000) {
        case 0b00010000:
            _compute_block_g<static_cast<Request>(0b00010000)>(block_begin, block_end, buffer, g, dg, hm);
            break;
        case 0b00100000:
            _compute_block_g<static_cast<Request>(0b00100000)>(block_begin, block_end, buffer, g, dg, hm);
            break;
        case 0b00110000:
            _compute_block_g<static_cast<Request>(0b00110000)>(block_begin, block_end, buffer, g, dg, hm);
            break;
        case 0b01000000:
            _compute_block_g<static_cast<Request>(0b01000000)>(block_begin, block_end, buffer, g, dg, hm);
            break;
        case 0b01010000:
            _compute_block_g<static_cast<Request>(0b01010000)>(block_begin, block_end, buffer, g, dg, hm);
            break;
        case 0b01100000:
            _compute_block_g<static_cast<Request>(0b01100000)>(block_begin, block_end, buffer, g, dg, hm);
            break;
        case 0b01110000:
            _compute_block_g<static_cast<Request>(0b01110000)>(block_begin, block_end, buffer, g, dg, hm);
            break;
        }
    }

    void _compute_g_async(std::mutex& mutex, Request request, std::atomic<index>& begin, const index end, Ref<Vector> g, Ref<Vector> dg, Ref<Vector> hm)
    {
        Vector future_buffer;
        Vector future_g;
        Vector future_dg;
        Vector future_hm;

        if (request & Request::G)
            future_g = Vector::Zero(len(g));

        if (request & Request::DG)
            future_dg = Vector::Zero(len(dg));

        if (request & Request::HG)
            future_hm = Vector::Zero(len(hm));

        if ((request & Request::G) || (request & Request::DG) || (request & Request::HG))
            future_buffer.resize(buffer_size());

        while (true) {
            index block_begin = begin.fetch_add(m_parallel_block_size);
            index block_end = std::min(block_begin + m_parallel_block_size, end);

            if (block_begin >= end)
                break;

            _compute_block_g(request, block_begin, block_end, future_buffer, future_g, future_dg, future_hm);
        }

        mutex.lock();

        if (request & Request::G)
            g += future_g;

        if (request & Request::DG)
            dg += future_dg;

        if (request & Request::HG)
            hm += future_hm;

        mutex.unlock();
    }

    // compute

    void _compute(Request request, Ref<Vector> buffer, double& f, Ref<Vector> g, Ref<Vector> df, Ref<Vector> dg, Ref<Vector> hm)
    {
        if (parallel_tasks() < 2) {
            _compute_block_f(request, 0, nb_objectives(), buffer, f, df, hm);

            if (sigma() != 1.0) {
                f *= m_sigma;
                df *= m_sigma;
                hm *= m_sigma;
            }

            _compute_block_g(request, 0, nb_constraints(), buffer, g, dg, hm);
        } else {
            std::vector<std::future<void>> tasks(parallel_tasks());

            std::mutex mutex;

            { // objectives
                std::atomic<index> begin = 0;

                for (index i = 0; i < len(tasks); i++)
                    tasks[i] = std::async(std::launch::async, &Problem::_compute_f_async, this, std::ref(mutex), request, std::ref(begin), nb_objectives(), std::ref(f), Ref<Vector>(df), Ref<Vector>(hm));

                for (index i = 0; i < len(tasks); i++)
                    tasks[i].wait();
            }

            if (sigma() != 1.0) {
                f *= m_sigma;
                df *= m_sigma;
                hm *= m_sigma;
            }

            { // constraints
                std::atomic<index> begin = 0;

                for (index i = 0; i < len(tasks); i++)
                    tasks[i] = std::async(std::launch::async, &Problem::_compute_g_async, this, std::ref(mutex), request, std::ref(begin), nb_constraints(), std::ref(g), Ref<Vector>(dg), Ref<Vector>(hm));

                for (index i = 0; i < len(tasks); i++)
                    tasks[i].wait();
            }
        }
    }

    double compute(Ref<Vector> g, Ref<Vector> df, Ref<Vector> dg, Ref<Vector> hm)
    {
        Request request = Request::F;

        if (len(g) != 0)
            request = static_cast<Request>(request | Request::G);

        if (len(df) != 0)
            request = static_cast<Request>(request | Request::DF);

        if (len(dg) != 0)
            request = static_cast<Request>(request | Request::DG);

        if (len(hm) != 0)
            request = static_cast<Request>(request | Request::HM);

        double f = 0.0;
        g.setZero();
        df.setZero();
        dg.setZero();
        hm.setZero();

        Ref<Vector> b = buffer();

        _compute(request, b, f, g, df, dg, hm);

        return f;
    }

    double compute(const Request request)
    {
        // double dummy;

        // Map<Vector> g(&dummy, 0);
        // Map<Vector> df(&dummy, 0);
        // Map<Vector> dg(&dummy, 0);
        // Map<Vector> hm(&dummy, 0);

        // if (request & Request::G) {
        //     g = this->g();
        // }

        Vector dummy(0);

        Ref<Vector> g(dummy);
        Ref<Vector> df(dummy);
        Ref<Vector> dg(dummy);
        Ref<Vector> hm(dummy);
        
        if (request & Request::G) {
            g = this->g();
        }

        if (request & Request::DF) {
            df = this->df();
        }

        if (request & Request::DG) {
            dg = dg_values();
        }
        
        if (request & Request::HM) {
            hm = hm_values();
        }

        return compute(g, df, dg, hm);
    }

    void eval(Ref<Vector> x)
    {
        set_x(x);

        f() = compute(g(), df(), dg_values(), hm_values());
    }

    // void eval_df(Ref<Vector> x)
    // {
    //     set_x(x);

    //     f() = compute(g(), df(), dg_values(), hm_values());
    // }

    // buffer_size

    index m_buffer_size;

    index buffer_size() const
    {
        return m_buffer_size;
    }

    // buffer

    Vector m_buffer;

    Ref<Vector> buffer()
    {
        if (len(m_buffer) != buffer_size())
            m_buffer.resize(buffer_size());

        return m_buffer;
    }

    // nb_nonzeros_dg

    index m_nb_nonzeros_dg;

    index nb_nonzeros_dg() const
    {
        return m_nb_nonzeros_dg;
    }

    // nb_nonzeros_hm

    index m_nb_nonzeros_hm;

    index nb_nonzeros_hm() const
    {
        return m_nb_nonzeros_hm;
    }

    // f

    double m_f;

    double& f()
    {
        return m_f;
    }

    // g

    Vector m_g;

    Ref<Vector> g()
    {
        if (len(m_g) != nb_equations())
            m_g.resize(nb_equations());

        return m_g;
    }

    // df

    Vector m_df;

    Ref<Vector> df()
    {
        if (len(m_df) != nb_variables())
            m_df.resize(nb_variables());

        return m_df;
    }

    // dg_values

    Vector m_dg_values;

    Ref<Vector> dg_values()
    {
        if (len(m_dg_values) != nb_nonzeros_dg())
            m_dg_values.resize(nb_nonzeros_dg());

        return m_dg_values;
    }

    // hm_values

    Vector m_hm_values;

    Ref<Vector> hm_values()
    {
        if (len(m_hm_values) != nb_nonzeros_hm())
            m_hm_values.resize(nb_nonzeros_hm());

        return m_hm_values;
    }

    // dg

    Map<const Sparse> dg() noexcept
    {
        return Map<const Sparse>(m_structure_dg.rows(), m_structure_dg.cols(), m_structure_dg.nb_nonzeros(), m_structure_dg.ia().data(), m_structure_dg.ja().data(), dg_values().data());
    }

    // hm

    Map<const Sparse> hm() noexcept
    {
        return Map<const Sparse>(m_structure_hm.rows(), m_structure_hm.cols(), m_structure_hm.nb_nonzeros(), m_structure_hm.ia().data(), m_structure_hm.ja().data(), hm_values().data());
    }

    // x

    Vector x() const
    {
        Vector result(nb_variables());

        for (index i = 0; i < nb_variables(); i++) {
            result(i) = m_variables[i]->value();
        }

        return result;
    }

    void set_x(const Ref<const Vector>& value)
    {
        if (len(value) != nb_variables()) {
            throw std::invalid_argument("invalid size");
        }

        for (index i = 0; i < nb_variables(); i++) {
            m_variables[i]->set_value(value(i));
        }
    }

    void add_x(const Ref<const Vector>& value)
    {
        if (len(value) != nb_variables()) {
            throw std::invalid_argument("invalid size");
        }

        for (index i = 0; i < nb_variables(); i++) {
            m_variables[i]->m_value += value(i);
        }
    }

    void sub_x(const Ref<const Vector>& value)
    {
        if (len(value) != nb_variables()) {
            throw std::invalid_argument("invalid size");
        }

        for (index i = 0; i < nb_variables(); i++) {
            m_variables[i]->m_value -= value(i);
        }
    }

    // x_lower_bounds

    Vector x_lower_bounds() const
    {
        Vector result(nb_variables());

        for (index i = 0; i < nb_variables(); i++)
            result(i) = m_variables[i]->lower_bound();

        return result;
    }

    void set_x_lower_bounds(const Ref<const Vector>& value)
    {
        if (len(value) != nb_variables())
            throw std::invalid_argument("invalid size");

        for (index i = 0; i < nb_variables(); i++)
            m_variables[i]->set_lower_bound(value(i));
    }

    // x_upper_bounds

    Vector x_upper_bounds() const
    {
        Vector result(nb_variables());

        for (index i = 0; i < nb_variables(); i++)
            result(i) = m_variables[i]->upper_bound();

        return result;
    }

    void set_x_upper_bounds(const Ref<const Vector>& value)
    {
        if (len(value) != nb_variables())
            throw std::invalid_argument("invalid size");

        for (index i = 0; i < nb_variables(); i++)
            m_variables[i]->set_upper_bound(value(i));
    }

    // parallel_tasks

    int m_parallel_tasks;

    int parallel_tasks() const
    {
        return m_parallel_tasks;
    }

    void set_parallel_tasks(const int value)
    {
        m_parallel_tasks = value;
    }

    // parallel_block_size

    index m_parallel_block_size;

    index parallel_block_size() const
    {
        return m_parallel_block_size;
    }

    void set_parallel_block_size(const index value)
    {
        m_parallel_block_size = value;
    }

    // linear_solver

    Pointer<LinearSolver> m_linear_solver;

    Pointer<LinearSolver> linear_solver()
    {
        if (m_linear_solver == nullptr) {
#ifdef EQLIB_USE_MKL
            m_linear_solver = new_<linear_solvers::PardisoLDLT>();
#else
            m_linear_solver = new_<linear_solvers::SimplicialLDLT>();
#endif
        }

        return m_linear_solver;
    }

    void set_linear_solver(const Pointer<LinearSolver> value)
    {
        m_linear_solver = value;
    }

    // df_norm
    
    double df_norm()
    {
        return m_df.norm();
    }

    // hm_inv_v
    
    Vector hm_inv_v(Ref<const Vector> v)
    {
        if (nb_variables() == 0) {
            return Vector(0);
        }

        if (linear_solver()->factorize(m_structure_hm.ia(), m_structure_hm.ja(), hm_values())) {
            throw std::runtime_error("Factorization failed");
        }

        Vector x(nb_variables());

        if (linear_solver()->solve(m_structure_hm.ia(), m_structure_hm.ja(), hm_values(), v, x)) {
            throw std::runtime_error("Solve failed");
        }

        return x;
    }

    // hm_v

    Vector hm_v(Ref<const Vector> v)
    {
        return hm().selfadjointView<Eigen::Upper>() * v.transpose();
    }

    // hm_diagonal

    Vector hm_diagonal()
    {
        Ref<Vector> hm = hm_values();

        Vector result(nb_variables());

        for (int row = 0; row < nb_variables(); row++) {
            const int i = m_structure_hm.ia(row);
            result(row) = hm(i);
        }

        return result;
    }

    void set_hm_diagonal(Eigen::Ref<const Vector> value)
    {
        Ref<Vector> hm = hm_values();

        for (int row = 0; row < nb_variables(); row++) {
            const int i = m_structure_hm.ia(row);
            hm(i) = value(row);
        }
    }

    void hm_add_diagonal(const double value)
    {
        Ref<Vector> hm = hm_values();

        for (int row = 0; row < nb_variables(); row++) {
            const int i = m_structure_hm.ia(row);
            hm(i) += value;
        }
    }

    // hm_norm_inf

    double hm_norm_inf()
    {
        Ref<Vector> hm = hm_values();

        Vector row_sum = Vector::Zero(nb_variables());

        for (int row = 0; row < nb_variables(); row++) {
            for (int i = m_structure_hm.ia(row); i < m_structure_hm.ia(row + 1); i++) {
                const int col = m_structure_hm.ja(i);

                const double abs_value = std::abs(hm(i));

                row_sum(row) += abs_value;

                if (row != col) {
                    row_sum(col) += abs_value;
                }
            }
        }

        return row_sum.maxCoeff();
    }

    // scale_hm

    void scale_hm(const double factor)
    {
        hm_values() *= factor;
    }

    // newton_step
    
    void newton_step()
    {
        if (nb_variables() == 0) {
            return;
        }

        if (linear_solver()->factorize(m_structure_hm.ia(), m_structure_hm.ja(), hm_values())) {
            throw std::runtime_error("Factorization failed");
        }

        Vector x(nb_variables());

        if (linear_solver()->solve(m_structure_hm.ia(), m_structure_hm.ja(), hm_values(), df(), x)) {
            throw std::runtime_error("Solve failed");
        }

        for (index i = 0; i < nb_variables(); i++) {
            auto& variable = *m_variables[i];
            
            variable.set_value(variable.value() - x[i]);
        }
    }
};

} // namespace eqlib
