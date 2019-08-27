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

private:    // variables
    Objectives m_objectives;
    Constraints m_constraints;

    double m_sigma;

    std::vector<Pointer<Equation>> m_equations;
    std::vector<Pointer<Variable>> m_variables;

    google::dense_hash_map<Pointer<Equation>, int> m_equation_indices;
    google::dense_hash_map<Pointer<Variable>, int> m_variable_indices;

    std::vector<std::tuple<std::vector<int>, std::vector<int>>> m_element_indices;

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
        tsl::robin_set<Pointer<Equation>> equation_set;
        tsl::robin_set<Pointer<Variable>> variable_set;

        // get lists of unique variables and equations

        for (auto it = m_objectives.begin(); it != m_objectives.end(); ++it) {
            for (const auto& variable : (*it)->variables()) {
                const auto find = variable_set.find(variable);
                
                if (find != variable_set.end()) {
                    continue;
                }

                variable_set.insert(variable);
                m_variables.push_back(variable);
            }
        }

        for (auto it = m_constraints.begin(); it != m_constraints.end(); ++it) {
            for (const auto& equation : (*it)->equations()) {
                const auto find = equation_set.find(equation);
                
                if (find != equation_set.end()) {
                    continue;
                }

                equation_set.insert(equation);
                m_equations.push_back(equation);
            }

            for (const auto& variable : (*it)->variables()) {
                const auto find = variable_set.find(variable);
                
                if (find != variable_set.end()) {
                    continue;
                }

                variable_set.insert(variable);
                m_variables.push_back(variable);
            }
        }

        // compute indices for variables and equations

        m_equation_indices.set_empty_key(nullptr);
        m_variable_indices.set_empty_key(nullptr);

        m_equation_indices.resize(m_equations.size());
        m_variable_indices.resize(m_variables.size());

        for (int i = 0; i < m_equations.size(); i++) {
            const auto& equation = m_equations[i];
            m_equation_indices[equation] = i;
        }

        for (int i = 0; i < m_variables.size(); i++) {
            const auto& variable = m_variables[i];
            m_variable_indices[variable] = i;
        }

        // compute indices for elements

        m_element_indices.reserve(m_objectives.size() + m_constraints.size());

        for (const auto& objective : m_objectives) {
            const auto variables = objective->variables();

            std::vector<int> variable_indices;
            variable_indices.reserve(variables.size());

            for (const auto& variable : variables) {
                variable_indices.push_back(m_variable_indices[variable]);
            }

            m_element_indices.emplace_back(std::vector<int>(0), variable_indices);
        }

        for (const auto& constraint : m_constraints) {
            const auto equations = constraint->equations();
            const auto variables = constraint->variables();

            std::vector<int> equation_indices;
            equation_indices.reserve(equations.size());

            std::vector<int> variable_indices;
            variable_indices.reserve(variables.size());

            for (const auto& equation : equations) {
                equation_indices.push_back(m_equation_indices[equation]);
            }

            for (const auto& variable : variables) {
                variable_indices.push_back(m_variable_indices[variable]);
            }

            m_element_indices.emplace_back(equation_indices, variable_indices);
        }

        // analyse sparse patterns

        const int n = m_variables.size();
        const int m = m_equations.size();

        std::vector<tsl::robin_set<int>> m_pattern_dg(n);
        std::vector<tsl::robin_set<int>> m_pattern_hl(n);

        for (int i = 0; i < m_objectives.size(); i++) {
            const auto& [equation_indices, variable_indices] = m_element_indices[i];

            for (const auto row : variable_indices) {
                for (const auto col : variable_indices) {
                    if (col > row) { // only lower triangle
                        continue;
                    }
                    m_pattern_hl[col].insert(row);
                }
            }
        }

        for (int i = 0; i < m_constraints.size(); i++) {
            const auto& [equation_indices, variable_indices] =
                m_element_indices[m_objectives.size() + i];

            for (const auto row : equation_indices) {
                for (const auto col : variable_indices) {
                    m_pattern_dg[col].insert(row);
                }
            }

            for (const auto row : variable_indices) {
                for (const auto col : variable_indices) {
                    if (col > row) { // only lower triangle
                        continue;
                    }
                    m_pattern_hl[col].insert(row);
                }
            }
        }

        Eigen::VectorXi size_dg(n);
        Eigen::VectorXi size_hl(n);

        for (int col = 0; col < n; col++) {
            size_dg[col] = m_pattern_dg[col].size();
            size_hl[col] = m_pattern_hl[col].size();
        }

        // allocate memory

        m_g = Vector(m);

        m_df = Vector(n);

        m_dg = Sparse(m, n);
        m_dg.reserve(size_dg);
        
        m_hl = Sparse(n, n);
        m_hl.reserve(size_hl);

        for (int col = 0; col < n; col++) {
            for (const int row : m_pattern_dg[col]) {
                m_dg.insert(row, col) = 1;
            }

            for (const int row : m_pattern_hl[col]) {
                m_hl.insert(row, col) = 1;
            }
        }
    }

public:     // methods
    void compute()
    {
        m_f = 0.0;
        m_g.setZero();
        m_df.setZero();
        Map<Vector>(m_dg.valuePtr(), m_dg.nonZeros()).setZero();
        Map<Vector>(m_hl.valuePtr(), m_hl.nonZeros()).setZero();
        
        for (int i = 0; i < m_objectives.size(); i++) {
            const auto& [equation_indices, variable_indices] =
                m_element_indices[i];

            const auto& objective = m_objectives[i];

            const auto n = variable_indices.size();
            
            Vector g(n);
            Matrix h(n, n);

            const double f = objective->compute(g, h);

            m_f += f;

            for (int local_row = 0; local_row < n; local_row++) {
                const int global_row = variable_indices[local_row];
                m_df[global_row] += g(local_row);

                for (int local_col = 0; local_col < n; local_col++) {
                    const int global_col = variable_indices[local_col];

                    if (global_col > global_row) {
                        continue;
                    }

                    m_hl.coeffRef(global_row, global_col) +=
                        h(local_row, local_col);
                }
            }
        }

        m_f *= sigma();
        m_df *= sigma();
        Map<Vector>(m_hl.valuePtr(), m_hl.nonZeros()) *= sigma();

        for (int i = 0; i < m_constraints.size(); i++) {
            const auto& [equation_indices, variable_indices] = m_element_indices[m_objectives.size() + i];

            const auto& constraint = m_constraints[i];

            const auto m = equation_indices.size();
            const auto n = variable_indices.size();

            std::vector<double> buffer(m * n + m * n * n);
            
            Vector rs(m);
            std::vector<Ref<Vector>> gs;
            std::vector<Ref<Matrix>> hs;
            
            gs.reserve(m);
            hs.reserve(m);
            
            for (int k = 0; k < m; k++) {
                Map<Vector> g(buffer.data() + k * n, n);
                Map<Matrix> h(buffer.data() + m * n + k * n * n, n, n);
                gs.push_back(g);
                hs.push_back(h);
            }

            constraint->compute(rs, gs, hs);

            for (int local_row = 0; local_row < m; local_row++) {
                const int global_row = equation_indices[local_row];
                m_g[global_row] += rs(local_row);

                for (int local_col = 0; local_col < n; local_col++) {
                    const int global_col = variable_indices[local_col];

                    m_dg.coeffRef(global_row, global_col) +=
                        gs[local_row](local_col);
                }
            }

            for (int k = 0; k < m; k++) {
                const auto& equation = m_equations[equation_indices[k]];

                hs[k] *= equation->multiplier();

                for (int local_row = 0; local_row < n; local_row++) {
                    const int global_row = variable_indices[local_row];

                    for (int local_col = 0; local_col < n; local_col++) {
                        const int global_col = variable_indices[local_col];

                        if (global_col > global_row) {
                            continue;
                        }

                        m_hl.coeffRef(global_row, global_col) +=
                            hs[k](local_row, local_col);
                    }
                }
            }
        }
    }

    Vector x() const
    {
        Vector result(nb_variables());

        for (int i = 0; i < result.size(); i++) {
            result(i) = variable(i)->act_value();
        }

        return result;
    }

    void set_x(Ref<const Vector> value) const
    {
        if (value.size() != nb_variables()) {
            throw std::runtime_error("Invalid size");
        }

        for (int i = 0; i < value.size(); i++) {
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

        for (int i = 0; i < result.size(); i++) {
            result(i) = variable(i)->multiplier();
        }

        return result;
    }

    void set_variable_multipliers(Ref<const Vector> value) const
    {
        if (value.size() != nb_variables()) {
            throw std::runtime_error("Invalid size");
        }

        for (int i = 0; i < value.size(); i++) {
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

        for (int i = 0; i < result.size(); i++) {
            result(i) = equation(i)->multiplier();
        }

        return result;
    }

    void set_equation_multipliers(Ref<const Vector> value) const
    {
        if (value.size() != nb_equations()) {
            throw std::runtime_error("Invalid size");
        }

        for (int i = 0; i < value.size(); i++) {
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

    int nb_equations() const noexcept
    {
        return static_cast<int>(m_equations.size());
    }

    int nb_variables() const noexcept
    {
        return static_cast<int>(m_variables.size());
    }

    const Pointer<Variable>& variable(const int index) const noexcept
    {
        return m_variables.at(index);
    }

    const Pointer<Equation>& equation(const int index) const noexcept
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

    double g(const int index) const noexcept
    {
        return m_g(index);
    }

    const Vector& df() const noexcept
    {
        return m_df;
    }

    double df(const int index) const noexcept
    {
        return m_df(index);
    }

    const Sparse& dg() const noexcept
    {
        return m_dg;
    }

    const Sparse& hl() const noexcept
    {
        return m_hl;
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
            .def_property_readonly("dg", &Type::dg)
            .def_property_readonly("hl", &Type::hl)
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