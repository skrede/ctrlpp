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
#include <cstdint>
#include <algorithm>

namespace ctrlpp::detail
{

// NLopt calls back into this object by address, so it is neither copied nor moved; its owner
// re-partitions a fresh one after a copy or move.
template <typename Scalar>
class nlopt_constraint_rows
{
public:
    nlopt_constraint_rows() : ineq_upper_indices_{}, ineq_lower_indices_{}, eq_indices_{}, problem_{nullptr}, values_{}, step_x_{}, minus_{}, plus_{} {}

    nlopt_constraint_rows(const nlopt_constraint_rows&) = delete;
    nlopt_constraint_rows& operator=(const nlopt_constraint_rows&) = delete;

    auto partition(const nlp_problem<Scalar>& problem, bool rejects_equalities) -> ctrlpp::expected<void, nlopt_setup_error>
    {
        problem_ = &problem;
        eq_indices_.clear();
        ineq_upper_indices_.clear();
        ineq_lower_indices_.clear();
        for(std::int32_t i = 0; i < problem.n_constraints; ++i)
        {
            if(problem.c_lower[i] == problem.c_upper[i])
                eq_indices_.push_back(i);
            else
            {
                if(std::isfinite(static_cast<double>(problem.c_upper[i])))
                    ineq_upper_indices_.push_back(i);
                if(std::isfinite(static_cast<double>(problem.c_lower[i])))
                    ineq_lower_indices_.push_back(i);
            }
        }
        if(rejects_equalities && !eq_indices_.empty())
            return ctrlpp::unexpected(nlopt_setup_error::incompatible_equality_constraints);
        return {};
    }

    void register_callbacks(nlopt::opt& opt, double tolerance)
    {
        if(!eq_indices_.empty())
            opt.add_equality_mconstraint(group_callback<&nlopt_constraint_rows::eq_indices_, &nlp_problem<Scalar>::c_lower, above>, this,
                                         std::vector<double>(eq_indices_.size(), tolerance));
        if(!ineq_upper_indices_.empty())
            opt.add_inequality_mconstraint(group_callback<&nlopt_constraint_rows::ineq_upper_indices_, &nlp_problem<Scalar>::c_upper, above>, this,
                                           std::vector<double>(ineq_upper_indices_.size(), tolerance));
        if(!ineq_lower_indices_.empty())
            opt.add_inequality_mconstraint(group_callback<&nlopt_constraint_rows::ineq_lower_indices_, &nlp_problem<Scalar>::c_lower, below>, this,
                                           std::vector<double>(ineq_lower_indices_.size(), tolerance));
    }

    auto max_violation(const std::vector<double>& x) const -> Scalar
    {
        if(problem_ == nullptr || problem_->n_constraints == 0)
            return Scalar{0};
        evaluate({x.data(), x.size()});
        Scalar max_viol{0};
        for(std::int32_t i = 0; i < problem_->n_constraints; ++i)
        {
            const auto ci = static_cast<Scalar>(values_[static_cast<std::size_t>(i)]);
            Scalar viol{0};
            if(ci > problem_->c_upper[i])
                viol = ci - problem_->c_upper[i];
            else if(ci < problem_->c_lower[i])
                viol = problem_->c_lower[i] - ci;
            max_viol = std::max(max_viol, viol);
        }
        return max_viol;
    }

private:
    std::vector<std::int32_t> ineq_upper_indices_;
    std::vector<std::int32_t> ineq_lower_indices_;
    std::vector<std::int32_t> eq_indices_;
    const nlp_problem<Scalar>* problem_;
    mutable std::vector<double> values_;
    mutable std::vector<double> step_x_;
    mutable std::vector<double> minus_;
    mutable std::vector<double> plus_;

    static auto above(double value, double bound) -> double
    {
        return value - bound;
    }

    static auto below(double value, double bound) -> double
    {
        return bound - value;
    }

    template <auto Indices, auto Bound, auto Residual>
    static void group_callback(unsigned m, double* result, unsigned n, const double* x, double* grad, void* data)
    {
        auto& self = *static_cast<nlopt_constraint_rows*>(data);
        const std::vector<std::int32_t>& indices = self.*Indices;
        const Eigen::VectorX<Scalar>& bound = self.problem_->*Bound;
        self.evaluate({x, n});
        for(unsigned i = 0; i < m; ++i)
        {
            const auto idx = indices[i];
            result[i] = Residual(self.values_[static_cast<std::size_t>(idx)], static_cast<double>(bound[idx]));
        }
        if(grad)
            self.template central_difference_jacobian<Residual>(m, n, x, grad, indices, bound);
    }

    template <auto Residual>
    void central_difference_jacobian(unsigned m, unsigned n, const double* x, double* grad, const std::vector<std::int32_t>& indices, const Eigen::VectorX<Scalar>& bound) const
    {
        const auto step_scale = std::cbrt(std::numeric_limits<double>::epsilon());
        step_x_.assign(x, x + n);
        plus_.resize(static_cast<std::size_t>(problem_->n_constraints));
        minus_.resize(static_cast<std::size_t>(problem_->n_constraints));
        for(unsigned j = 0; j < n; ++j)
        {
            const double h_raw = step_scale * std::max(1.0, std::abs(x[j]));
            const double temp = x[j] + h_raw;
            const double h = temp - x[j];
            evaluate_shifted(n, j, step_x_[j] + h, plus_);
            evaluate_shifted(n, j, step_x_[j] - h, minus_);
            for(unsigned i = 0; i < m; ++i)
            {
                const auto idx = static_cast<std::size_t>(indices[i]);
                const double row_bound = static_cast<double>(bound[indices[i]]);
                grad[static_cast<std::size_t>(i) * n + j] = (Residual(plus_[idx], row_bound) - Residual(minus_[idx], row_bound)) / (2.0 * h);
            }
        }
    }

    void evaluate_shifted(unsigned n, unsigned j, double shifted, std::vector<double>& out) const
    {
        const double original = step_x_[j];
        step_x_[j] = shifted;
        problem_->constraints_cb(std::span<const double>{step_x_.data(), n}, std::span<double>{out.data(), out.size()});
        step_x_[j] = original;
    }

    void evaluate(std::span<const double> x) const
    {
        values_.resize(static_cast<std::size_t>(problem_->n_constraints));
        problem_->constraints_cb(x, std::span<double>{values_.data(), values_.size()});
    }
};

}

#endif
