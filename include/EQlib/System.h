#pragma once

#include "Assemble.h"
#include "Define.h"
#include "Dof.h"
#include "Element.h"
#include "Log.h"
#include "Timer.h"

#include <Eigen/PardisoSupport>

#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace EQlib {

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
    std::vector<Dof> m_dofs;

    std::unordered_map<Dof, int> m_dof_indices;

    int m_nb_free_dofs;
    int m_nb_fixed_dofs;

    std::vector<std::shared_ptr<Element>> m_elements;
    std::vector<std::vector<Index>> m_index_table;

    std::vector<std::pair<std::shared_ptr<Element>, std::vector<Index>>> m_element_index_table;

    std::vector<std::vector<int>> m_pattern;
    Eigen::VectorXi m_col_nonzeros;

    Sparse m_lhs;
    Vector m_rhs;
    Vector m_x;
    Vector m_target;
    Vector m_residual;

    int m_stopping_reason;

    Eigen::PardisoLDLT<Sparse, Eigen::Upper> m_solver;

    double m_load_factor;

    int m_nb_threads;

public:     // constructors
    System(
        std::vector<std::shared_ptr<Element>> elements)
    : m_load_factor(1)
    , m_nb_threads(0)
    {
        initialize(std::move(elements));
    }

private:    // methods
    void initialize(
        std::vector<std::shared_ptr<Element>> elements)
    {
        Log log;
        log.info(1, "==> Initialize system...");
        log.info(2, "The system consists of {} elements", elements.size());

        Timer timer;

        // get dofs

        log.info(3, "Getting element dofs...");

        const auto nb_elements = elements.size();

        std::vector<std::vector<Dof>> element_dofs(nb_elements);

        for (size_t i = 0; i < nb_elements; i++) {
            const auto& element = elements[i];

            element_dofs[i] = element->dofs();
        }

        // create set of unique dofs

        log.info(3, "Creating set of unique dofs...");

        std::unordered_set<Dof> dof_set;
        std::vector<Dof> free_dofs;
        std::vector<Dof> fixed_dofs;

        for (size_t i = 0; i < nb_elements; i++) {
            const auto& element = elements[i];

            for (const auto& dof : element_dofs[i]) {
                if (dof_set.find(dof) != dof_set.end()) {
                    continue;
                }

                dof_set.insert(dof);

                if (dof.isfixed()) {
                    fixed_dofs.push_back(dof);
                } else {
                    free_dofs.push_back(dof);
                }
            }
        }

        // store system size

        m_nb_free_dofs = static_cast<int>(free_dofs.size());
        m_nb_fixed_dofs = static_cast<int>(fixed_dofs.size());

        log.info(2, "The system has {} free and {} fixed dofs", m_nb_free_dofs,
            m_nb_fixed_dofs);

        // create a vector of unique dofs

        m_dofs.reserve(m_nb_free_dofs + m_nb_fixed_dofs);

        m_dofs.insert(m_dofs.end(), free_dofs.begin(), free_dofs.end());
        m_dofs.insert(m_dofs.end(), fixed_dofs.begin(), fixed_dofs.end());

        // create a {dof -> index} map

        log.info(3, "Creating dof indices...");

        for (int i = 0; i < m_dofs.size(); i++) {
            m_dof_indices[m_dofs[i]] = i;
        }

        // create index table

        for (size_t i = 0; i < nb_elements; i++) {
            const auto& element = elements[i];

            const auto& dofs = element_dofs[i];

            const auto nb_dofs = dofs.size();

            std::vector<Index> dof_indices(nb_dofs);

            for (int local = 0; local < nb_dofs; local++) {
                const auto dof = dofs[local];
                const auto global = m_dof_indices[dof];
                dof_indices[local] = {local, global};
            }

            std::sort(dof_indices.begin(), dof_indices.end());

            m_index_table.push_back(dof_indices);
            m_element_index_table.push_back(std::make_pair(element, dof_indices));
        }

        // analyze pattern

        log.info(3, "Analyzing pattern...");

        std::vector<std::unordered_set<int>> pattern(m_nb_free_dofs);

        for (const auto& dof_indices : m_index_table) {
            const size_t nb_dofs = dof_indices.size();

            for (size_t row = 0; row < nb_dofs; row++) {
                const auto row_index = dof_indices[row];

                if (row_index.global >= nb_free_dofs()) {
                    continue;
                }

                for (size_t col = row; col < nb_dofs; col++) {
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

        m_elements = std::move(elements);

        // setup system vectors and matrices

        log.info(3, "Allocating memory...");

        m_lhs = Sparse(nb_free_dofs(), nb_free_dofs());

        if (nb_free_dofs() > 0) {
            m_lhs.reserve(m_col_nonzeros);

            for (int col = 0; col < m_pattern.size(); col++) {
                for (const int row : m_pattern[col]) {
                    m_lhs.insert(row, col);
                }
            }
        }

        log.info(2, "The system matrix has {} nonzero entries ({:.3f}%)", m_lhs.nonZeros(), m_lhs.nonZeros() * 100.0 / m_lhs.size());

        m_rhs = Vector(nb_free_dofs());

        m_x = Vector(nb_free_dofs());

        m_target = Vector(nb_free_dofs());
        m_residual = Vector(nb_free_dofs());

        // setup solver

        m_solver.analyzePattern(m_lhs);

        log.info(1, "System initialized in {:.3f} sec", timer.ellapsed());
    }

public:     // getters and setters
    const std::vector<Dof>& dofs() const
    {
        return m_dofs;
    }

    int nb_dofs() const
    {
        return static_cast<int>(m_dofs.size());
    }

    int nb_free_dofs() const
    {
        return m_nb_free_dofs;
    }

    int nb_fixed_dofs() const
    {
        return m_nb_fixed_dofs;
    }

    Sparse lhs() const
    {
        return m_lhs;
    }

    Vector rhs() const
    {
        return m_rhs;
    }

    Vector x() const
    {
        Vector result(nb_free_dofs());

        for (int i = 0; i < result.size(); i++) {
            result[i] = m_dofs[i].delta();
        }

        return result;
    }

    void set_x(Ref<const Vector> value) const
    {
        if (value.size() != nb_free_dofs()) {
            throw std::runtime_error("Invalid size");
        }

        for (int i = 0; i < value.size(); i++) {
            m_dofs[i].set_delta(value[i]);
        }
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

    int nb_threads() const
    {
        return m_nb_threads;
    }

    void set_nb_threads(const int value)
    {
        m_nb_threads = value;
    }

    std::vector<std::shared_ptr<Element>> elements() const
    {
        return m_elements;
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
    int dof_index(const Dof& dof) const
    {
        return m_dof_indices.at(dof);
    }

    void compute()
    {
        Log log;
        log.info(1, "==> Computing system...");

        Timer timer;

        // compute lhs and rhs

        Assemble::parallel(m_nb_threads, m_element_index_table, m_lhs, m_rhs);

        log.info(1, "System computed in {:.3f} sec", timer.ellapsed());
    }

    void solve()
    {
        // options

        const int maxiter = 100;
        const double rtol = 1e-7;
        const double xtol = 1e-7;
        const bool parallel = true;

        // setup

        Log log;
        log.info(1, "==> Solving system...");

        Timer timer;

        for (int i = 0; i < nb_free_dofs(); i++) {
            m_target[i] = m_dofs[i].target();
        }


        if (m_load_factor != 1.0) {
            m_target *= m_load_factor;
        }

        for (int iteration = 0; ; iteration++) {
            // check max iterations

            if (iteration >= maxiter) {
                m_stopping_reason = 2;
                log.info(2, "Stopped because iteration >= {}", maxiter);
                break;
            }

            log.info(2, "Iteration {}", iteration);

            // compute lhs and rhs

            log.info(2, "Computing system...");

            Assemble::parallel(m_nb_threads, m_element_index_table, m_lhs, m_rhs);

            // check residual

            log.info(2, "Computing residual...");

            m_residual = m_target + m_rhs;

            const double rnorm = m_residual.norm();

            log.info(2, "The norm of the residual is {}", rnorm);

            // check residual norm

            if (rnorm < rtol) {
                m_stopping_reason = 0;
                log.info(2, "Stopped because rnorm < {}", rtol);
                break;
            }

            // solve iteration

            log.info(2, "Solving the linear equation system...");

            m_solver.factorize(m_lhs);

            if (!m_solver.info() == Eigen::Success) {
                throw std::runtime_error("Factorization failed");
            }

            m_x = m_solver.solve(m_residual);

            if (!m_solver.info() == Eigen::Success) {
                throw std::runtime_error("Solve failed");
            }

            // update system

            log.info(2, "Updating system...");

            for (int i = 0; i < nb_free_dofs(); i++) {
                m_dofs[i].set_delta(m_dofs[i].delta() + m_x(i));
            }

            // check x norm

            const double xnorm = m_x.norm();

            log.info(2, "The norm of the step is {}", xnorm);

            if (xnorm < xtol) {
                log.info(2, "Stopped because xnorm < {}", xtol);
                m_stopping_reason = 1;
                break;
            }
        }

        for (int i = 0; i < nb_free_dofs(); i++) {
            m_dofs[i].set_residual(m_residual(i));
        }

        log.info(1, "System solved in {:.3f} sec", timer.ellapsed());
    }

    void solve_linear()
    {
        // setup

        for (int i = 0; i < nb_dofs(); i++) {
            m_target[i] = m_dofs[i].target();
        }

        if (m_load_factor != 1.0) {
            m_target *= m_load_factor;
        }

        m_residual = m_target + m_rhs;

        // compute lhs and rhs

        compute();

        // solve

        m_solver.factorize(m_lhs);

        if (!m_solver.info() == Eigen::Success) {
            throw std::runtime_error("Factorization failed");
        }

        m_x = m_solver.solve(m_residual);

        if (!m_solver.info() == Eigen::Success) {
            throw std::runtime_error("Solve failed");
        }

        // update system

        for (int i = 0; i < nb_free_dofs(); i++) {
            m_dofs[i].set_delta(m_x(i));
        }
    }
};

} // namespace EQlib