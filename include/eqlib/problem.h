#pragma once

#include "common.h"
#include "objective.h"
#include "constraint.h"
#include "variable.h"
#include "sparse_structure.h"

#include <algorithm>
#include <future>
#include <stdexcept>

namespace eqlib
{
    struct Index
    {
        index local;
        index global;

        bool operator<(const Index& other) const noexcept
        {
            return global < other.global;
        }
    };

    struct ObjectiveIndices
    {
        ObjectiveIndices() = default;

        ObjectiveIndices(const Pointer<Objective>& objective, const index nb_variables, const index nb_constants) : m_objective(objective), variable_indices(nb_variables), constant_indices(nb_constants)
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

    struct ConstraintIndices
    {
        ConstraintIndices() = default;

        ConstraintIndices(const Pointer<Constraint>& constraint, const index nb_equations, const index nb_variables, const index nb_constants) : m_constraint(constraint), equation_indices(nb_equations), variable_indices(nb_variables), constant_indices(nb_constants)
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

    struct Problem
    {
        Problem() : m_buffer_size(0), m_nb_nonzeros_dg(0), m_nb_nonzeros_hm(0), m_parallel_level(0), m_parallel_block_size(256)
        {
        }

        static Problem create(const std::vector<Pointer<Objective>>& objectives, const std::vector<Pointer<Constraint>>& constraints)
        {
            Logger::info("Create a new system");

            Problem result;

            const index nb_objectives = len(objectives);
            const index nb_constraints = len(constraints);

            Logger::info("{} objective elements", nb_objectives);
            Logger::info("{} constraints elements", nb_constraints);

            result.m_objective_indices = std::vector<ObjectiveIndices>(nb_objectives);
            result.m_constraint_indices = std::vector<ConstraintIndices>(nb_constraints);

            // index variables and equations

            for (index i = 0; i < nb_objectives; i++)
            {
                const auto& objective = objectives[i];

                index local_nb_constants = 0;

                for (const auto& variable : objective->variables())
                {
                    if (variable->is_active())
                    {
                        result.m_variable_indices.try_emplace(variable, len(result.m_variable_indices));
                    }
                    else
                    {
                        local_nb_constants += 1;
                        result.m_constant_indices.try_emplace(variable, len(result.m_constant_indices));
                    }
                }

                const index local_nb_variables = objective->nb_variables() - local_nb_constants;

                result.m_objective_indices[i] = ObjectiveIndices(objective, local_nb_variables, local_nb_constants);
            }

            for (index i = 0; i < nb_constraints; i++)
            {
                const auto& constraint = constraints[i];

                index local_nb_equations = 0;

                for (const auto& equation : constraint->equations())
                {
                    if (equation->is_active())
                    {
                        local_nb_equations += 1;
                        result.m_equation_indices.try_emplace(equation, len(result.m_equation_indices));
                    }
                }

                index local_nb_constants = 0;

                for (const auto& variable : constraint->variables())
                {
                    if (variable->is_active())
                    {
                        result.m_variable_indices.try_emplace(variable, len(result.m_variable_indices));
                    }
                    else
                    {
                        local_nb_constants += 1;
                        result.m_constant_indices.try_emplace(variable, len(result.m_constant_indices));
                    }
                }

                const index local_nb_variables = constraint->nb_variables() - local_nb_constants;

                result.m_constraint_indices[i] = ConstraintIndices(constraint, local_nb_equations, local_nb_variables, local_nb_constants);
            }

            // collect equations

            index nb_equations = len(result.m_equation_indices);

            Logger::info("{} equations", nb_equations);

            result.m_equations = Equations(nb_equations);

            for (const auto& [equation, index] : result.m_equation_indices)
            {
                result.m_equations[index] = equation;
            }

            // collect variables

            index nb_variables = len(result.m_variable_indices);

            Logger::info("{} variables", nb_variables);

            result.m_variables = Variables(nb_variables);

            for (const auto& [variable, index] : result.m_variable_indices)
            {
                result.m_variables[index] = variable;
            }

            // collect constants

            index nb_constants = len(result.m_constant_indices);

            result.m_constants = Variables(nb_constants);

            for (const auto& [constant, index] : result.m_constant_indices)
            {
                result.m_constants[index] = constant;
            }

            // assign indices

            for (auto& objective_indices : result.m_objective_indices)
            {
                const auto& objective = objective_indices.objective();

                auto variable_indices_it = objective_indices.variable_indices.begin();
                auto constant_indices_it = objective_indices.constant_indices.begin();

                for (index local_index = 0; local_index < objective.nb_variables(); local_index++)
                {
                    const auto variable = objective.variable(local_index);

                    if (variable->is_active())
                    {
                        const index global_index = result.m_variable_indices[variable];
                        *(variable_indices_it++) = { local_index, global_index };
                    }
                    else
                    {
                        const index global_index = result.m_constant_indices[variable];
                        *(constant_indices_it++) = { local_index, global_index };
                    }
                }

                objective_indices.sort();
            }

            for (auto& constraint_indices : result.m_constraint_indices)
            {
                const auto& constraint = constraint_indices.constraint();

                auto equation_indices_it = constraint_indices.equation_indices.begin();

                for (index local_index = 0; local_index < constraint.nb_equations(); local_index++)
                {
                    const auto equation = constraint.equation(local_index);

                    if (equation->is_active())
                    {
                        const index global_index = result.m_equation_indices[equation];
                        *(equation_indices_it++) = { local_index, global_index };
                    }
                }

                auto variable_indices_it = constraint_indices.variable_indices.begin();
                auto constant_indices_it = constraint_indices.constant_indices.begin();

                for (index local_index = 0; local_index < constraint.nb_variables(); local_index++)
                {
                    const auto variable = constraint.variable(local_index);

                    if (variable->is_active())
                    {
                        const index global_index = result.m_variable_indices[variable];
                        *(variable_indices_it++) = { local_index, global_index };
                    }
                    else
                    {
                        const index global_index = result.m_constant_indices[variable];
                        *(constant_indices_it++) = { local_index, global_index };
                    }
                }

                constraint_indices.sort();
            }

            // pattern

            std::vector<RobinSet<index>> pattern_dg(nb_equations);
            std::vector<RobinSet<index>> pattern_hm(nb_variables);

            for (auto& constraint_indices : result.m_constraint_indices)
            {
                const auto& equation_indices = constraint_indices.equation_indices;
                const auto& variable_indices = constraint_indices.variable_indices;

                for (auto& index_i : equation_indices)
                {
                    for (auto& index_j : variable_indices)
                    {
                        pattern_dg[index_i.global].insert(index_j.global);
                    }
                }

                const index nb_variables = len(variable_indices);

                for (index i = 0; i < nb_variables; i++)
                {
                    const index global_i = variable_indices[i].global;

                    for (index j = i; j < nb_variables; j++)
                    {
                        const index global_j = variable_indices[j].global;

                        pattern_hm[global_i].insert(global_j);
                    }
                }

                // update buffer size

                const auto& constraint = constraint_indices.constraint();

                const index required_buffer_size = constraint.nb_equations() + constraint.nb_equations() * constraint.nb_variables() + constraint.nb_variables() * constraint.nb_variables();

                if (required_buffer_size > result.m_buffer_size)
                {
                    result.m_buffer_size = required_buffer_size;
                }
            }

            for (auto& objective_indices : result.m_objective_indices)
            {
                const index nb_variables = len(objective_indices.variable_indices);

                // update pattern

                for (index i = 0; i < nb_variables; i++)
                {
                    const index global_i = objective_indices.variable_indices[i].global;

                    for (index j = i; j < nb_variables; j++)
                    {
                        const index global_j = objective_indices.variable_indices[j].global;

                        pattern_hm[global_i].insert(global_j);
                    }
                }

                // update buffer size

                const auto& objective = objective_indices.objective();

                const index required_buffer_size = objective.nb_variables() + objective.nb_variables() * objective.nb_variables();

                if (required_buffer_size > result.m_buffer_size)
                {
                    result.m_buffer_size = required_buffer_size;
                }
            }

            for (const auto& i : pattern_dg)
            {
                result.m_nb_nonzeros_dg += len(i);
            }

            for (const auto& i : pattern_hm)
            {
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

            // index variables and equations

            for (index i = 0; i < nb_objectives; i++)
            {
                const auto& objective = objectives[i];

                index local_nb_constants = 0;

                for (const auto& variable : objective->variables())
                {
                    if (!variable->is_active())
                    {
                        local_nb_constants += 1;
                    }
                }

                const index local_nb_variables = objective->nb_variables() - local_nb_constants;

                result.m_objective_indices[i] = ObjectiveIndices(objective, local_nb_variables, local_nb_constants);
            }

            for (index i = 0; i < nb_constraints; i++)
            {
                const auto& constraint = constraints[i];

                index local_nb_equations = 0;

                for (const auto& equation : constraint->equations())
                {
                    if (equation->is_active())
                    {
                        local_nb_equations += 1;
                    }
                }

                index local_nb_constants = 0;

                for (const auto& variable : constraint->variables())
                {
                    if (!variable->is_active())
                    {
                        local_nb_constants += 1;
                    }
                }

                const index local_nb_variables = constraint->nb_variables() - local_nb_constants;

                result.m_constraint_indices[i] = ConstraintIndices(constraint, local_nb_equations, local_nb_variables, local_nb_constants);
            }

            // assign indices

            for (auto& objective_indices : result.m_objective_indices)
            {
                const auto& objective = objective_indices.objective();

                auto variable_indices_it = objective_indices.variable_indices.begin();
                auto constant_indices_it = objective_indices.constant_indices.begin();

                for (index local_index = 0; local_index < objective.nb_variables(); local_index++)
                {
                    const auto variable = objective.variable(local_index);

                    if (variable->is_active())
                    {
                        const auto it = prototype.m_variable_indices.find(variable);

                        if (it != prototype.m_variable_indices.end())
                        {
                            const index global_index = it.value();
                            *(variable_indices_it++) = { local_index, global_index };
                        }
                    }
                    else
                    {
                        const auto it = prototype.m_constant_indices.find(variable);

                        if (it != prototype.m_constant_indices.end())
                        {
                            const index global_index = it.value();
                            *(constant_indices_it++) = { local_index, global_index };
                        }
                    }
                }

                objective_indices.sort();
            }

            for (auto& constraint_indices : result.m_constraint_indices)
            {
                const auto& constraint = constraint_indices.constraint();

                auto equation_indices_it = constraint_indices.equation_indices.begin();

                for (index local_index = 0; local_index < constraint.nb_equations(); local_index++)
                {
                    const auto equation = constraint.equation(local_index);

                    if (equation->is_active())
                    {
                        const auto it = prototype.m_equation_indices.find(equation);

                        if (it != prototype.m_equation_indices.end())
                        {
                            const index global_index = it.value();
                            *(equation_indices_it++) = { local_index, global_index };
                        }
                    }
                }

                auto variable_indices_it = constraint_indices.variable_indices.begin();
                auto constant_indices_it = constraint_indices.constant_indices.begin();

                for (index local_index = 0; local_index < constraint.nb_variables(); local_index++)
                {
                    const auto variable = constraint.variable(local_index);

                    if (variable->is_active())
                    {
                        const auto it = prototype.m_variable_indices.find(variable);

                        if (it != prototype.m_variable_indices.end())
                        {
                            const index global_index = it.value();
                            *(variable_indices_it++) = { local_index, global_index };
                        }
                    }
                    else
                    {
                        const auto it = prototype.m_constant_indices.find(variable);

                        if (it != prototype.m_constant_indices.end())
                        {
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

        RobinMap<Pointer<Variable>, index> m_variable_indices;

        index variable_index(const Pointer<Variable>& variable) const
        {
            const auto it = m_variable_indices.find(variable);

            if (it == m_variable_indices.end())
                return -1;

            return it.value();
        }

        // constant indices

        RobinMap<Pointer<Variable>, index> m_constant_indices;

        index constant_index(const Pointer<Variable>& variable) const
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

        Variables m_variables;

        index nb_variables() const
        {
            return len(m_variables);
        }

        const Variables& variables() const
        {
            return m_variables;
        }

        // constants

        Variables m_constants;

        index nb_constants() const
        {
            return len(m_constants);
        }

        const Variables& constants() const
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

        // compute f

        template <Request R>
        void _compute_block_f(const index block_begin, const index block_end, Ref<Vector> buffer, double& f, Ref<Vector> df, Ref<Vector> hm)
        {
            constexpr bool request_f = (R & Request::F) != 0;
            constexpr bool request_df = (R & Request::Df) != 0;
            constexpr bool request_hf = (R & Request::Hf) != 0;

            for (index i = block_begin; i < block_end; i++)
            {
                auto& objective_indices = m_objective_indices[i];

                auto& objective = objective_indices.objective();

                if (!objective.is_active())
                    return;

                const index nb_variables = objective.nb_variables();

                Map<Vector> local_df(buffer.data(), nb_variables);
                Map<Matrix> local_hm(buffer.data() + nb_variables, nb_variables, nb_variables);

                const double local_f = objective.compute(local_df, local_hm, R);

                if constexpr (request_f)
                    f += local_f;

                if constexpr (request_df)
                {
                    for (index i = 0; i < len(objective_indices.variable_indices); i++)
                    {
                        const auto index_i = objective_indices.variable_indices[i];

                        df(index_i.global) += local_df(index_i.local);
                    }
                }

                if constexpr (request_hf)
                {
                    for (index i = 0; i < len(objective_indices.variable_indices); i++)
                    {
                        const auto index_i = objective_indices.variable_indices[i];

                        for (index j = i; j < len(objective_indices.variable_indices); j++)
                        {
                            const auto index_j = objective_indices.variable_indices[j];

                            const index idx = m_structure_hm.get_index(index_i.global, index_j.global);

                            assert(idx >= 0);

                            hm(idx) += local_hm(index_i.local, index_j.local);
                        }
                    }
                }
            }
        }

        template <Request R>
        void _compute_f(std::atomic<index>& begin, const index end, Ref<Vector> buffer, double& f, Ref<Vector> df, Ref<Vector> hm, const index level)
        {
            if (end - begin <= 0)
                return;

            constexpr bool request_f = (R & Request::F) != 0;
            constexpr bool request_df = (R & Request::Df) != 0;
            constexpr bool request_hf = (R & Request::Hf) != 0;

            if constexpr (!request_f && !request_df && !request_hf)
                return;

            if (level < m_parallel_level)
            {
                Vector future_buffer = Vector::Zero(len(buffer));
                double future_f = 0.0;
                Vector future_df = Vector::Zero(len(df));
                Vector future_hm = Vector::Zero(len(hm));

                auto future = std::async(std::launch::async, &Problem::_compute_f<R>, this, std::ref(begin), end, Ref<Vector>(future_buffer), std::ref(future_f), Ref<Vector>(future_df), Ref<Vector>(future_hm), level + 1);

                _compute_f<R>(begin, end, buffer, f, df, hm, level + 1);

                future.wait();

                f += future_f;
                df += future_df;
                hm += future_hm;

                return;
            }

            if (m_parallel_level == 0) // serial
            {
                _compute_block_f<R>(begin, end, buffer, f, df, hm);
            }
            else // parallel
            {
                while (true)
                {
                    const index block_begin = begin.fetch_add(m_parallel_block_size);

                    if (block_begin >= end)
                        break;

                    const index block_end = std::min(block_begin + m_parallel_block_size, end);

                    _compute_block_f<R>(block_begin, block_end, buffer, f, df, hm);
                }
            }
        }

        template <Request R>
        void compute_f(Ref<Vector> buffer, double& f, Ref<Vector> df, Ref<Vector> hm)
        {
            std::atomic<index> i;

            _compute_f<R>(i, len(m_objective_indices), buffer, f, df, hm, 0);
        }

        void compute_f(Ref<Vector> buffer, double& f, Ref<Vector> df, Ref<Vector> hm)
        {
            if (len(df) == 0 && len(hm) == 0)
            {
                constexpr Request R = static_cast<Request>(Request::F);
                compute_f<R>(buffer, f, EmptyVector, EmptyVector);
            }
            else if (len(hm) == 0)
            {
                constexpr Request R = static_cast<Request>(Request::F | Request::Df);
                compute_f<R>(buffer, f, df, EmptyVector);
            }
            else if (len(df) == 0)
            {
                constexpr Request R = static_cast<Request>(Request::F | Request::Hf);
                compute_f<R>(buffer, f, EmptyVector, hm);
            }
            else
            {
                constexpr Request R = static_cast<Request>(Request::F | Request::Df | Request::Hf);
                compute_f<R>(buffer, f, df, hm);
            }
        }

        void compute_f(double& f, Ref<Vector> df, Ref<Vector> hm)
        {
            Ref<Vector> b = buffer();
            compute_f(b, f, df, hm);
        }

        // compute g

        template <Request R>
        void _compute_block_g(const index block_begin, const index block_end, Ref<Vector> buffer, Ref<Vector> g, Ref<Vector> dg, Ref<Vector> hm)
        {
            constexpr bool request_g = (R & Request::G) != 0;
            constexpr bool request_dg = (R & Request::Dg) != 0;
            constexpr bool request_hm = (R & Request::Hg) != 0;

            if constexpr (!request_g && !request_dg && !request_hm)
                return;

            for (index i = block_begin; i < block_end; i++)
            {
                auto& constraint_indices = m_constraint_indices[i];

                auto& constraint = constraint_indices.constraint();

                if (!constraint.is_active())
                    continue;

                const index nb_equations = constraint.nb_equations();
                const index nb_variables = constraint.nb_variables();

                Map<Vector> local_g(buffer.data(), nb_equations);
                Map<Matrix> local_dg(buffer.data() + nb_equations, nb_equations, nb_variables);
                Map<Matrix> local_hm(buffer.data() + nb_equations + nb_equations * nb_variables, nb_variables, nb_variables);

                constraint.compute(local_g, local_dg, local_hm, R);

                if constexpr (request_g)
                {
                    for (index i = 0; i < len(constraint_indices.equation_indices); i++)
                    {
                        const auto index_i = constraint_indices.equation_indices[i];

                        if constexpr (request_g)
                            g(index_i.global) += local_g(index_i.local);

                        if constexpr (request_dg)
                        {
                            for (index j = 0; j < len(constraint_indices.variable_indices); j++)
                            {
                                const auto index_j = constraint_indices.variable_indices[j];

                                const index idx = m_structure_dg.get_index(index_i.global, index_j.global);

                                assert(idx >= 0);

                                dg(idx) += local_dg(index_i.local, index_j.local);
                            }
                        }
                    }

                    if constexpr (request_hm)
                    {
                        for (index i = 0; i < len(constraint_indices.variable_indices); i++)
                        {
                            const auto index_i = constraint_indices.variable_indices[i];

                            for (index j = i; j < len(constraint_indices.variable_indices); j++)
                            {
                                const auto index_j = constraint_indices.variable_indices[j];

                                const index idx = m_structure_hm.get_index(index_i.global, index_j.global);

                                assert(idx >= 0);

                                hm(idx) += local_hm(index_i.local, index_j.local);
                            }
                        }
                    }
                }
            }
        }

        template <Request R>
        void _compute_g(std::atomic<index>& begin, const index end, Ref<Vector> buffer, Ref<Vector> g, Ref<Vector> dg, Ref<Vector> hm, const index level)
        {
            if (end - begin <= 0)
                return;

            constexpr bool request_g = (R & Request::G) != 0;
            constexpr bool request_dg = (R & Request::Dg) != 0;
            constexpr bool request_hm = (R & Request::Hg) != 0;

            if constexpr (!request_g && !request_dg && !request_hm)
                return;

            if (level < m_parallel_level)
            {
                Vector future_buffer = Vector::Zero(len(buffer));
                Vector future_g = Vector::Zero(len(g));
                Vector future_dg = Vector::Zero(len(dg));
                Vector future_hm = Vector::Zero(len(hm));

                auto future = std::async(std::launch::async, &Problem::_compute_g<R>, this, std::ref(begin), end, Ref<Vector>(future_buffer), Ref<Vector>(future_g), Ref<Vector>(future_dg), Ref<Vector>(future_hm), level + 1);

                _compute_g<R>(begin, end, buffer, g, dg, hm, level + 1);

                future.wait();

                g += future_g;
                dg += future_dg;
                hm += future_hm;

                return;
            }

            if (m_parallel_level == 0) // serial
            {
                _compute_block_g<R>(begin, end, buffer, g, dg, hm);
            }
            else // parallel
            {
                while (true)
                {
                    const index block_begin = begin.fetch_add(m_parallel_block_size);

                    if (block_begin >= end)
                        break;

                    const index block_end = std::min(block_begin + m_parallel_block_size, end);

                    _compute_block_g<R>(block_begin, block_end, buffer, g, dg, hm);
                }
            }
        }

        template <Request R>
        void compute_g(Ref<Vector> buffer, Ref<Vector> g, Ref<Vector> dg, Ref<Vector> hm)
        {
            std::atomic<index> i;

            _compute_g<R>(i, len(m_constraint_indices), buffer, g, dg, hm, 0);
        }

        void compute_g(Ref<Vector> buffer, Ref<Vector> g, Ref<Vector> dg, Ref<Vector> hm)
        {
            if (len(g) == 0 && len(dg) == 0 && len(hm) == 0)
            {
                return;
            }
            else if (len(dg) == 0 && len(hm) == 0)
            {
                constexpr Request R = static_cast<Request>(Request::G);
                compute_g<R>(buffer, g, EmptyVector, EmptyVector);
            }
            else if (len(g) == 0 && len(hm) == 0)
            {
                constexpr Request R = static_cast<Request>(Request::Dg);
                compute_g<R>(buffer, EmptyVector, dg, EmptyVector);
            }
            else if (len(g) == 0 && len(dg) == 0)
            {
                constexpr Request R = static_cast<Request>(Request::Hg);
                compute_g<R>(buffer, EmptyVector, EmptyVector, hm);
            }
            else if (len(hm) == 0)
            {
                constexpr Request R = static_cast<Request>(Request::G | Request::Dg);
                compute_g<R>(buffer, g, dg, EmptyVector);
            }
            else if (len(dg) == 0)
            {
                constexpr Request R = static_cast<Request>(Request::G | Request::Hg);
                compute_g<R>(buffer, g, EmptyVector, hm);
            }
            else if (len(g) == 0)
            {
                constexpr Request R = static_cast<Request>(Request::Dg | Request::Hg);
                compute_g<R>(buffer, EmptyVector, dg, hm);
            }
            else
            {
                constexpr Request R = static_cast<Request>(Request::G | Request::Dg | Request::Hg);
                compute_g<R>(buffer, g, dg, hm);
            }
        }

        void compute_g(Ref<Vector> g, Ref<Vector> dg, Ref<Vector> hm)
        {
            Ref<Vector> b = buffer();
            compute_g(b, g, dg, hm);
        }

        // compute

        template <Request R>
        double compute(Ref<Vector> g, Ref<Vector> df, Ref<Vector> dg, Ref<Vector> hm)
        {
            double f = 0.0;
            g.setZero();
            df.setZero();
            dg.setZero();
            hm.setZero();

            Ref<Vector> b = buffer();

            compute_f<R>(b, f, df, hm);
            compute_g<R>(b, g, dg, hm);

            return f;
        }

        double compute(Ref<Vector> g, Ref<Vector> df, Ref<Vector> dg, Ref<Vector> hm)
        {
            double f = 0.0;
            g.setZero();
            df.setZero();
            dg.setZero();
            hm.setZero();

            Ref<Vector> b = buffer();

            compute_f(b, f, df, hm);
            compute_g(b, g, dg, hm);

            return f;
        }

        void eval(Ref<Vector> x)
        {
            set_x(x);

            f() = compute(g(), df(), dg(), hm());
        }

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

        // dg

        Vector m_dg;

        Ref<Vector> dg()
        {
            if (len(m_dg) != nb_nonzeros_dg())
                m_dg.resize(nb_nonzeros_dg());

            return m_dg;
        }

        // hm

        Vector m_hm;

        Ref<Vector> hm()
        {
            if (len(m_hm) != nb_nonzeros_hm())
                m_hm.resize(nb_nonzeros_hm());

            return m_hm;
        }

        // x

        Vector x() const
        {
            Vector result(nb_variables());

            for (index i = 0; i < nb_variables(); i++)
            {
                result(i) = m_variables[i]->value();
            }

            return result;
        }

        void set_x(const Ref<const Vector>& value)
        {
            if (len(value) != nb_variables())
            {
                throw std::invalid_argument("invalid size");
            }

            for (index i = 0; i < nb_variables(); i++)
            {
                m_variables[i]->set_value(value(i));
            }
        }

        void add_x(const Ref<const Vector>& value)
        {
            if (len(value) != nb_variables())
            {
                throw std::invalid_argument("invalid size");
            }

            for (index i = 0; i < nb_variables(); i++)
            {
                m_variables[i]->m_value += value(i);
            }
        }

        void sub_x(const Ref<const Vector>& value)
        {
            if (len(value) != nb_variables())
            {
                throw std::invalid_argument("invalid size");
            }

            for (index i = 0; i < nb_variables(); i++)
            {
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

        // parallel_level

        index m_parallel_level;

        index parallel_level() const
        {
            return m_parallel_level;
        }

        void set_parallel_level(const index value)
        {
            m_parallel_level = value;
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
    };

} // namespace eqlib
