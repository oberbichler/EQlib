#pragma once

#include "Define.h"
#include "Constraint.h"
#include "Objective.h"
#include "Settings.h"

#include <sparsehash/dense_hash_map>

#include <tsl/robin_set.h>

#include <tuple>
#include <utility>
#include <vector>

namespace EQlib {

class Problem
{
private:    // types
    using Objectives = std::vector<Pointer<Objective>>;
    using Constraints = std::vector<Pointer<Constraint>>;

    using Equations = std::vector<Pointer<Equation>>;
    using Variables = std::vector<Pointer<Variable>>;

    struct Index
    {
        Index(index local, index global)
        : local(local), global(global)
        {
        }

        index local;
        index global;

        bool operator<(const Index& other) const noexcept
        {
            return global < other.global;
        }
    };

private:    // variables
    Objectives m_objectives;
    Constraints m_constraints;

    double m_sigma;

    std::vector<Pointer<Equation>> m_equations;
    std::vector<Pointer<Variable>> m_variables;

    google::dense_hash_map<Pointer<Equation>, index> m_equation_indices;
    google::dense_hash_map<Pointer<Variable>, index> m_variable_indices;

    std::vector<index> m_element_f_nb_variables;
    std::vector<index> m_element_g_nb_variables;
    std::vector<index> m_element_g_nb_equations;

    std::vector<std::vector<Index>> m_element_f_variable_indices;
    std::vector<std::vector<Index>> m_element_f_df_indices;

    std::vector<std::vector<Index>> m_element_g_equation_indices;
    std::vector<std::vector<Index>> m_element_g_variable_indices;

    double m_f;
    Vector m_g;
    Vector m_df;
    Sparse m_dg;
    Sparse m_hl;

public:     // constructors
    Problem(
        Objectives objectives,
        Constraints constraints,
        Settings linear_solver)
    : m_objectives(std::move(objectives))
    , m_constraints(std::move(constraints))
    , m_sigma(1.0)
    {
        Log::info(1, "==> Initialize problem...");

        const auto nb_elements_f = length(m_objectives);
        const auto nb_elements_g = length(m_constraints);

        Log::info(2, "The problem consists of {} objective elements", nb_elements_f);
        Log::info(2, "The problem consists of {} constraint elements", nb_elements_g);


        Log::info(3, "Getting equations and variables...");

        m_element_f_nb_variables.resize(nb_elements_f);
        m_element_g_nb_variables.resize(nb_elements_g);
        m_element_g_nb_equations.resize(nb_elements_g);

        std::vector<Variables> variables_f(nb_elements_f);

        std::vector<Equations> equations_g(nb_elements_g);
        std::vector<Variables> variables_g(nb_elements_g);

        for (index i = 0; i < nb_elements_f; i++) {
            const auto& element = *m_objectives[i];

            variables_f[i] = element.variables();

            m_element_f_nb_variables[i] = length(variables_f[i]);
        }

        for (index i = 0; i < nb_elements_g; i++) {
            const auto& element = *m_constraints[i];

            equations_g[i] = element.equations();
            variables_g[i] = element.variables();

            m_element_g_nb_variables[i] = length(variables_g[i]);
            m_element_g_nb_equations[i] = length(equations_g[i]);
        }


        Log::info(3, "Creating the set of unique equations...");

        tsl::robin_set<Pointer<Equation>> equation_set;

        for (const auto& equations : equations_g) {
            for (const auto& equation : equations) {
                const auto find = equation_set.find(equation);

                if (find != equation_set.end()) {
                    continue;
                }

                equation_set.insert(equation);
                m_equations.push_back(equation);
            }
        }


        Log::info(3, "Creating the set of unique variables...");

        tsl::robin_set<Pointer<Variable>> variable_set;

        for (const auto& variables : variables_f) {
            for (const auto& variable : variables) {
                if (variable->is_fixed()) {
                    continue;
                }

                const auto find = variable_set.find(variable);

                if (find != variable_set.end()) {
                    continue;
                }

                variable_set.insert(variable);
                m_variables.push_back(variable);
            }
        }

        for (const auto& variables : variables_g) {
            for (const auto& variable : variables) {
                if (variable->is_fixed()) {
                    continue;
                }

                const auto find = variable_set.find(variable);

                if (find != variable_set.end()) {
                    continue;
                }

                variable_set.insert(variable);
                m_variables.push_back(variable);
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

                if (variable->is_fixed()) {
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

                if (variable->is_fixed()) {
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

        std::vector<tsl::robin_set<index>> m_pattern_dg(n);
        std::vector<tsl::robin_set<index>> m_pattern_hl(n);

        for (index i = 0; i < length(m_objectives); i++) {
            const auto& variable_indices = m_element_f_variable_indices[i];

            for (index col_i = 0; col_i < length(variable_indices); col_i++) {
                const auto col = variable_indices[col_i];

                for (index row_i = col_i; row_i < length(variable_indices); row_i++) {
                    const auto row = variable_indices[row_i];

                    m_pattern_hl[col.global].insert(row.global);
                }
            }
        }

        for (index i = 0; i < length(m_constraints); i++) {
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

        Eigen::VectorXi sparse_size_dg(n);
        Eigen::VectorXi sparse_size_hl(n);

        for (index col = 0; col < n; col++) {
            sparse_size_dg[col] = length(m_pattern_dg[col]);
            sparse_size_hl[col] = length(m_pattern_hl[col]);
        }


        Log::info(3, "Allocate memory...");

        m_g = Vector(m);

        m_df = Vector(n);

        if (m > 0 && n > 0) {
            m_dg = Sparse(m, n);
            m_dg.reserve(sparse_size_dg);

            for (index col = 0; col < n; col++) {
                for (const index row : m_pattern_dg[col]) {
                    m_dg.insert(row, col) = 1;
                }
            }
        }

        if (n > 0) {
            m_hl = Sparse(n, n);
            m_hl.reserve(sparse_size_hl);

            index i = 0;

            for (index col = 0; col < n; col++) {
                for (const index row : m_pattern_hl[col]) {
                    m_hl.insert(row, col) = 1;
                }
            }

        }
    }

public:     // methods
    void compute()
    {
        m_f = 0.0;
        m_g.setZero();
        m_df.setZero();
        dg_values().setZero();
        hl_values().setZero();

        for (index i = 0; i < length(m_objectives); i++) {
            const auto& variable_indices = m_element_f_variable_indices[i];

            const auto& objective = m_objectives[i];

            const auto n = m_element_f_nb_variables[i];

            Vector g(n);
            Matrix h(n, n);

            const double f = objective->compute(g, h);

            m_f += f;

            for (index col_i = 0; col_i < length(variable_indices); col_i++) {
                const auto col = variable_indices[col_i];

                df(col.global) += g(col.local);

                for (index row_i = col_i; row_i < length(variable_indices); row_i++) {
                    const auto row = variable_indices[row_i];

                    hl(row.global, col.global) += h(row.local, col.local);
                }
            }
        }

        m_f *= sigma();
        m_df *= sigma();
        hl_values() *= sigma();

        for (index i = 0; i < length(m_constraints); i++) {
            const auto& equation_indices = m_element_g_equation_indices[i];
            const auto& variable_indices = m_element_g_variable_indices[i];

            const auto& constraint = m_constraints[i];

            const auto m = m_element_g_nb_equations[i];
            const auto n = m_element_g_nb_variables[i];

            if (n == 0 || m == 0) {
                continue;
            }

            std::vector<double> buffer(m * n + m * n * n);

            Vector fs(m);
            std::vector<Ref<Vector>> gs;
            std::vector<Ref<Matrix>> hs;

            gs.reserve(m);
            hs.reserve(m);

            for (index k = 0; k < m; k++) {
                Map<Vector> g(buffer.data() + k * n, n);
                Map<Matrix> h(buffer.data() + m * n + k * n * n, n, n);
                gs.push_back(g);
                hs.push_back(h);
            }

            constraint->compute(fs, gs, hs);

            for (const auto& equation_index : equation_indices) {
                const auto& equation = m_equations[equation_index.global];

                g(equation_index.global) += fs(equation_index.local);

                auto& local_g = gs[equation_index.local];
                auto& local_h = hs[equation_index.local];

                local_h *= equation->multiplier();

                for (index col_i = 0; col_i < length(variable_indices); col_i++) {
                    const auto col = variable_indices[col_i];

                    dg(equation_index.global, col.global) += local_g(col.local);

                    for (index row_i = col_i; row_i < length(variable_indices); row_i++) {
                        const auto row = variable_indices[row_i];

                        hl(row.global, col.global) += local_h(row.local, col.local);
                    }
                }
            }
        }
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

    double f() const noexcept
    {
        return m_f;
    }

    const Vector& g() const noexcept
    {
        return m_g;
    }

    double& g(const index index)
    {
        return m_g(index);
    }

    const Vector& df() const noexcept
    {
        return m_df;
    }

    double& df(const index index)
    {
        return m_df(index);
    }

    const Sparse& dg() const noexcept
    {
        return m_dg;
    }

    Ref<Vector> dg_values() noexcept
    {
        return Map<Vector>(m_dg.valuePtr(), m_dg.nonZeros());
    }

    double& dg(const index index)
    {
        return *(m_dg.valuePtr() + index);
    }

    double& dg(const index row, const index col)
    {
        return m_dg.coeffRef(row, col);
    }

    const Sparse& hl() const noexcept
    {
        return m_hl;
    }

    Ref<Vector> hl_values() noexcept
    {
        return Map<Vector>(m_hl.valuePtr(), m_hl.nonZeros());
    }

    double& hl(const index index)
    {
        return *(m_hl.valuePtr() + index);
    }

    double& hl(const index row, const index col)
    {
        return m_hl.coeffRef(row, col);
    }

public:     // python
    template <typename TModule>
    static void register_python(TModule& m)
    {
        namespace py = pybind11;
        using namespace pybind11::literals;

        using Type = Problem;
        using Holder = Pointer<Type>;

        const std::string name = "Problem";

        py::class_<Type, Holder>(m, name.c_str())
            // constructors
            .def(py::init<Objectives, Constraints, Settings>(),
                "objectives"_a, "constraints"_a, "linear_solver"_a = Settings())
            // read-only properties
            .def_property_readonly("equations", &Type::equations)
            .def_property_readonly("variables", &Type::variables)
            .def_property_readonly("f", &Type::f)
            .def_property_readonly("g", py::overload_cast<>(&Type::g, py::const_))
            .def_property_readonly("df", py::overload_cast<>(&Type::df, py::const_))
            .def_property_readonly("dg", py::overload_cast<>(&Type::dg, py::const_))
            .def_property_readonly("hl", py::overload_cast<>(&Type::hl, py::const_))
            .def_property_readonly("x", &Type::x)
            .def_property_readonly("nb_equations", &Type::nb_equations)
            .def_property_readonly("nb_variables", &Type::nb_variables)
            .def_property_readonly("variable_multipliers", &Type::variable_multipliers)
            .def_property_readonly("equation_multipliers", &Type::equation_multipliers)
            // properties
            .def_property("sigma", &Type::sigma, &Type::set_sigma)
            // methods
            .def("compute", &Type::compute)
        ;
    }
};

} // namespace EQlib