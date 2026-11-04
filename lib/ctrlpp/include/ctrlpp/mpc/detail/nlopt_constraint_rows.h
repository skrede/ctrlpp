#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_NLOPT_CONSTRAINT_ROWS_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_NLOPT_CONSTRAINT_ROWS_H

#include "ctrlpp/expected.h"

#include "ctrlpp/mpc/nlp_types.h"
#include "ctrlpp/mpc/nlp_solver.h"

#include <nlopt.hpp>

#include <span>
#include <cmath>
#include <limits>
#include <vector>
#include <cstddef>
#include <algorithm>

namespace ctrlpp::detail
{

// NLopt calls back into this object by address, so it is neither copied nor moved; its owner
// re-partitions a fresh one after a copy or move.
template <typename Scalar>
class nlopt_constraint_rows
{
public:
    nlopt_constraint_rows()
        : m_ineq_upper_indices{}
        , m_ineq_lower_indices{}
        , m_eq_indices{}
        , m_problem{nullptr}
        , m_values{}
        , m_step_x{}
        , m_minus{}
        , m_plus{}
    {
    }

    nlopt_constraint_rows(const nlopt_constraint_rows&) = delete;
    nlopt_constraint_rows& operator=(const nlopt_constraint_rows&) = delete;

    ctrlpp::expected<void, nlopt_setup_error> partition(const nlp_problem<Scalar>& problem, bool rejects_equalities)
    {
        m_problem = &problem;
        m_eq_indices.clear();
        m_ineq_upper_indices.clear();
        m_ineq_lower_indices.clear();
        for(Eigen::Index i = 0; i < problem.n_constraints; ++i)
        {
            if(problem.c_lower[i] == problem.c_upper[i])
                m_eq_indices.push_back(i);
            else
            {
                if(std::isfinite(static_cast<double>(problem.c_upper[i])))
                    m_ineq_upper_indices.push_back(i);
                if(std::isfinite(static_cast<double>(problem.c_lower[i])))
                    m_ineq_lower_indices.push_back(i);
            }
        }
        if(rejects_equalities && !m_eq_indices.empty())
            return ctrlpp::unexpected(nlopt_setup_error::incompatible_equality_constraints);
        return {};
    }

    void register_callbacks(nlopt::opt& opt, double tolerance)
    {
        if(!m_eq_indices.empty())
            opt.add_equality_mconstraint(
                    group_callback<&nlopt_constraint_rows::m_eq_indices, &nlp_problem<Scalar>::c_lower, above>, this, std::vector<double>(m_eq_indices.size(), tolerance));
        if(!m_ineq_upper_indices.empty())
            opt.add_inequality_mconstraint(group_callback<&nlopt_constraint_rows::m_ineq_upper_indices, &nlp_problem<Scalar>::c_upper, above>, this,
                    std::vector<double>(m_ineq_upper_indices.size(), tolerance));
        if(!m_ineq_lower_indices.empty())
            opt.add_inequality_mconstraint(group_callback<&nlopt_constraint_rows::m_ineq_lower_indices, &nlp_problem<Scalar>::c_lower, below>, this,
                    std::vector<double>(m_ineq_lower_indices.size(), tolerance));
    }

    Scalar max_violation(const std::vector<double>& x) const
    {
        if(m_problem == nullptr || m_problem->n_constraints == 0)
            return Scalar{0};
        evaluate({x.data(), x.size()});
        Scalar max_viol{0};
        for(Eigen::Index i = 0; i < m_problem->n_constraints; ++i)
        {
            const auto ci = static_cast<Scalar>(m_values[static_cast<std::size_t>(i)]);
            Scalar viol{0};
            if(ci > m_problem->c_upper[i])
                viol = ci - m_problem->c_upper[i];
            else if(ci < m_problem->c_lower[i])
                viol = m_problem->c_lower[i] - ci;
            max_viol = std::max(max_viol, viol);
        }
        return max_viol;
    }

private:
    std::vector<Eigen::Index> m_ineq_upper_indices;
    std::vector<Eigen::Index> m_ineq_lower_indices;
    std::vector<Eigen::Index> m_eq_indices;
    // Null until partition(). An optional reference here grows the solve frame on Cortex-M7.
    const nlp_problem<Scalar>* m_problem;
    mutable std::vector<double> m_values;
    mutable std::vector<double> m_step_x;
    mutable std::vector<double> m_minus;
    mutable std::vector<double> m_plus;

    static double above(double value, double bound) { return value - bound; }

    static double below(double value, double bound) { return bound - value; }

    template <auto Indices, auto Bound, auto Residual>
    static void group_callback(unsigned m, double* result, unsigned n, const double* x, double* grad, void* data)
    {
        auto& self = *static_cast<nlopt_constraint_rows*>(data);
        const std::vector<Eigen::Index>& indices = self.*Indices;
        const Eigen::VectorX<Scalar>& bound = self.m_problem->*Bound;
        self.evaluate({x, n});
        for(unsigned i = 0; i < m; ++i)
        {
            const auto idx = indices[i];
            result[i] = Residual(self.m_values[static_cast<std::size_t>(idx)], static_cast<double>(bound[idx]));
        }
        if(grad)
            self.template central_difference_jacobian<Residual>(m, n, x, grad, indices, bound);
    }

    template <auto Residual>
    void central_difference_jacobian(unsigned m, unsigned n, const double* x, double* grad, const std::vector<Eigen::Index>& indices, const Eigen::VectorX<Scalar>& bound) const
    {
        const auto step_scale = std::cbrt(std::numeric_limits<double>::epsilon());
        m_step_x.assign(x, x + n);
        m_plus.resize(static_cast<std::size_t>(m_problem->n_constraints));
        m_minus.resize(static_cast<std::size_t>(m_problem->n_constraints));
        for(unsigned j = 0; j < n; ++j)
        {
            const double h_raw = step_scale * std::max(1.0, std::abs(x[j]));
            const double temp = x[j] + h_raw;
            const double h = temp - x[j];
            evaluate_shifted(n, j, m_step_x[j] + h, m_plus);
            evaluate_shifted(n, j, m_step_x[j] - h, m_minus);
            for(unsigned i = 0; i < m; ++i)
            {
                const auto idx = static_cast<std::size_t>(indices[i]);
                const double row_bound = static_cast<double>(bound[indices[i]]);
                grad[static_cast<std::size_t>(i) * n + j] = (Residual(m_plus[idx], row_bound) - Residual(m_minus[idx], row_bound)) / (2.0 * h);
            }
        }
    }

    void evaluate_shifted(unsigned n, unsigned j, double shifted, std::vector<double>& out) const
    {
        const double original = m_step_x[j];
        m_step_x[j] = shifted;
        m_problem->constraints_cb(std::span<const double>{m_step_x.data(), n}, std::span<double>{out.data(), out.size()});
        m_step_x[j] = original;
    }

    void evaluate(std::span<const double> x) const
    {
        m_values.resize(static_cast<std::size_t>(m_problem->n_constraints));
        m_problem->constraints_cb(x, std::span<double>{m_values.data(), m_values.size()});
    }
};

}

#endif
