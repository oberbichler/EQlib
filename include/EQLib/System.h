#pragma once

#include "Define.h"
#include "Element.h"
#include "Dof.h"

#include <fmt/format.h>

#include <unordered_map>
#include <unordered_set>
#include <set>

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

    std::vector<std::set<int>> m_pattern;
    Eigen::VectorXi m_row_nonzeros;

    Sparse m_lhs;
    Vector m_rhs;
    Vector m_x;

    int m_stopping_reason;

    SparseSolver m_solver;

public:     // methods
    System(
        std::vector<std::shared_ptr<Element>> elements,
        py::dict options)
    {
        // options

        const std::string linear_solver = get_or_default<std::string>(options,
            "linear_solver", "eigen_lu");
        const bool symmetric = get_or_default(options, "symmetric", true);

        // get dofs

        const auto nb_elements = elements.size();

        std::vector<std::vector<Dof>> element_dofs(nb_elements);

        for (size_t i = 0; i < nb_elements; i++) {
            const auto& element = elements[i];

            element_dofs[i] = element->dofs();
        }

        // create set of unique dofs

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

        // create a vector of unique dofs

        m_dofs.reserve(m_nb_free_dofs + m_nb_fixed_dofs);

        m_dofs.insert(m_dofs.end(), free_dofs.begin(), free_dofs.end());
        m_dofs.insert(m_dofs.end(), fixed_dofs.begin(), fixed_dofs.end());

        // create a {dof -> index} map

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
        }

        // analyze pattern

        std::vector<std::set<int>> pattern(m_nb_free_dofs);

        for (const auto& dof_indices : m_index_table) {
            for (const auto col_index : dof_indices) {
                if (col_index.global >= m_nb_free_dofs) {
                    continue;
                }

                auto hint = pattern[col_index.global].begin();

                const int max_col = symmetric ? col_index.global + 1 : m_nb_free_dofs;

                for (const auto row_index : dof_indices) {
                    if (row_index.global >= max_col) {
                        continue;
                    }

                    hint = pattern[col_index.global].insert(hint,
                        row_index.global);
                }
            }
        }

        m_row_nonzeros = Eigen::VectorXi(m_nb_free_dofs);

        for (int i = 0; i < pattern.size(); i++) {
            m_row_nonzeros[i] = static_cast<int>(pattern[i].size());
        }

        // store data

        m_pattern = std::move(pattern);

        m_elements = std::move(elements);

        // setup system vectors and matrices

        m_lhs = Sparse(nb_free_dofs(), nb_free_dofs());

        m_lhs.reserve(m_row_nonzeros);

        for (int col = 0; col < m_pattern.size(); col++) {
            for (const int row : m_pattern[col]) {
                m_lhs.insert(row, col);
            }
        }

        m_rhs = Vector(nb_free_dofs());

        m_x = Vector(nb_free_dofs());

        // setup solver

        m_solver.analyzePattern(m_lhs);
    }

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

    int dof_index(const Dof& dof) const
    {
        return m_dof_indices.at(dof);
    }

    Sparse lhs() const
    {
        return m_lhs;
    }

    Vector rhs() const
    {
        return m_rhs;
    }

    void compute(py::dict options)
    {
        // options

        const bool symmetric = get_or_default(options, "symmetric", true);

        // set lhs and rhs to zero

        for (int i = 0; i < m_lhs.outerSize(); i++) {
            for (Sparse::InnerIterator it(m_lhs, i); it; ++it){
                it.valueRef() = 0;
            }
        }

        m_rhs.setZero();

        // compute and add local lhs and rhs

        for (size_t i = 0; i < m_elements.size(); i++) {
            const auto& element = m_elements[i];

            const auto& dof_indices = m_index_table[i];

            const auto [local_lhs, local_rhs] = element->compute(options);

            for (const auto row_index : dof_indices) {
                if (row_index.global >= nb_free_dofs()) {
                    continue;
                }

                m_rhs(row_index.global) += local_rhs(row_index.local);

                const int max_col = symmetric ? row_index.global + 1 : nb_free_dofs();

                for (const auto col_index : dof_indices) {
                    if (col_index.global >= max_col) {
                        continue;
                    }

                    m_lhs.coeffRef(row_index.global, col_index.global) +=
                        local_lhs(row_index.local, col_index.local);
                }
            }
        }
    }

    void solve(py::dict options)
    {
        // options

        const double lambda = get_or_default(options, "lambda", 1.0);
        const int maxiter = get_or_default(options, "maxiter", 100);
        const double rtol = get_or_default(options, "rtol", 1e-7);
        const double xtol = get_or_default(options, "xtol", 1e-7);

        // setup

        Vector target(nb_free_dofs());
        Vector residual(nb_free_dofs());

        for (int i = 0; i < nb_free_dofs(); i++) {
            target[i] = m_dofs[i].target();
        }

        if (lambda != 1.0) {
            target *= lambda;
        }

        for (int iteration = 0; iteration < maxiter; iteration++) {
            // compute lhs and rhs

            compute(options);

            // check residual

            residual = rhs() - target;

            const double rnorm = residual.norm();

            if (rnorm < rtol) {
                std::cout << fmt::format("{:>4} {:}", iteration, rnorm) << std::endl;
                m_stopping_reason = 0;
                break;
            }

            //

            if (iteration + 1 == maxiter) {
                std::cout << fmt::format("{:>4} {:}", iteration, rnorm) << std::endl;
                m_stopping_reason = 2;
            } else {
                std::cout << fmt::format("{:>4} {:}", iteration, rnorm) << std::endl;
            }

            // solve iteration

            m_solver.factorize(m_lhs);
            m_x = m_solver.solve(m_rhs);

            const double xnorm = m_x.norm();

            if (xnorm < xtol) {
                std::cout << fmt::format("{:>4} {:}", iteration, rnorm) << std::endl;
                m_stopping_reason = 1;
                break;
            }

            for (int i = 0; i < nb_free_dofs(); i++) {
                m_dofs[i].set_delta(m_dofs[i].delta() - m_x(i));
            }
        }

        for (int i = 0; i < nb_free_dofs(); i++) {
            m_dofs[i].set_residual(residual(i));
        }
    }

    std::string stopping_reason_message() const
    {
        std::string message;

        switch (m_stopping_reason) {
        case -1:
            message = "Not solved";
            break;
        case 0:
            message = "A solution was found, given rtol";
            break;
        case 1:
            message = "A solution was found, given xtol";
            break;
        case 2:
            message = "The iteration limit was reached";
            break;
        default:
            message = "Error. Unknown stopping reason";
            break;
        }

        return message;
    }
};

} // namespace EQlib