#ifndef HPP_GUARD_CTRLPP_MPC_ARGMIN_PROBLEM_H
#define HPP_GUARD_CTRLPP_MPC_ARGMIN_PROBLEM_H

#include "ctrlpp/mpc/nlp_solver.h"

#include <Eigen/Core>

#include <span>
#include <cmath>
#include <vector>
#include <cstddef>
#include <algorithm>

namespace ctrlpp
{

template <typename Scalar>
class argmin_problem
{
public:
    static constexpr int problem_dimension = Eigen::Dynamic;

    void bind(const nlp_problem<Scalar>& prob)
    {
        problem_ = &prob;
    }

    auto value(const Eigen::VectorX<Scalar>& x) const -> Scalar
    {
        return problem_->cost(
            std::span<const Scalar>{x.data(), static_cast<std::size_t>(x.size())});
    }

    auto dimension() const -> int
    {
        return problem_->n_vars;
    }

    void gradient(const Eigen::VectorX<Scalar>& x, Eigen::VectorX<Scalar>& g) const
    {
        problem_->gradient(
            std::span<const Scalar>{x.data(), static_cast<std::size_t>(x.size())},
            std::span<Scalar>{g.data(), static_cast<std::size_t>(g.size())});
    }

    auto lower_bounds() const -> Eigen::VectorX<Scalar>
    {
        return problem_->x_lower;
    }

    auto upper_bounds() const -> Eigen::VectorX<Scalar>
    {
        return problem_->x_upper;
    }

protected:
    const nlp_problem<Scalar>* problem_{nullptr};
};

template <typename Scalar>
class argmin_constrained_problem : public argmin_problem<Scalar>
{
    using base = argmin_problem<Scalar>;

public:
    using base::problem_dimension;

    void partition(const nlp_problem<Scalar>& prob)
    {
        base::bind(prob);

        eq_indices.clear();
        ineq_upper_indices.clear();
        ineq_lower_indices.clear();

        for(int i = 0; i < prob.n_constraints; ++i)
        {
            if(prob.c_lower[i] == prob.c_upper[i])
            {
                eq_indices.push_back(i);
            }
            else
            {
                if(std::isfinite(static_cast<double>(prob.c_upper[i])))
                    ineq_upper_indices.push_back(i);
                if(std::isfinite(static_cast<double>(prob.c_lower[i])))
                    ineq_lower_indices.push_back(i);
            }
        }

        n_eq = static_cast<int>(eq_indices.size());
        n_ineq_upper = static_cast<int>(ineq_upper_indices.size());
        n_ineq_lower = static_cast<int>(ineq_lower_indices.size());

        raw_buf_.resize(static_cast<std::size_t>(prob.n_constraints));
        fd_x_buf_.resize(prob.n_vars);
    }

    void constraints(const Eigen::VectorX<Scalar>& x, Eigen::VectorX<Scalar>& c_out) const
    {
        eval_raw(x);

        int idx = 0;

        for(auto i : eq_indices)
            c_out[idx++] = raw_buf_[static_cast<std::size_t>(i)] - base::problem_->c_lower[i];

        for(auto i : ineq_upper_indices)
            c_out[idx++] = raw_buf_[static_cast<std::size_t>(i)] - base::problem_->c_upper[i];

        for(auto i : ineq_lower_indices)
            c_out[idx++] = base::problem_->c_lower[i] - raw_buf_[static_cast<std::size_t>(i)];
    }

    auto num_equality() const -> int
    {
        return n_eq;
    }

    auto num_inequality() const -> int
    {
        return n_ineq_upper + n_ineq_lower;
    }

    void constraint_jacobian(const Eigen::VectorX<Scalar>& x, Eigen::MatrixX<Scalar>& J) const
    {
        const int m = n_eq + n_ineq_upper + n_ineq_lower;
        const int n = base::problem_->n_vars;
        J.resize(m, n);

        if(base::problem_->constraint_jacobian)
        {
            Eigen::MatrixX<Scalar> J_raw(base::problem_->n_constraints, n);
            base::problem_->constraint_jacobian(
                std::span<const Scalar>{x.data(), static_cast<std::size_t>(x.size())},
                std::span<Scalar>{J_raw.data(), static_cast<std::size_t>(J_raw.size())});

            int idx = 0;
            for(auto i : eq_indices)
                J.row(idx++) = J_raw.row(i);
            for(auto i : ineq_upper_indices)
                J.row(idx++) = J_raw.row(i);
            for(auto i : ineq_lower_indices)
                J.row(idx++) = -J_raw.row(i);

            return;
        }

        const auto step_scale = std::cbrt(std::numeric_limits<Scalar>::epsilon());
        Eigen::VectorX<Scalar> c_plus(m);
        Eigen::VectorX<Scalar> c_minus(m);

        fd_x_buf_ = x;

        for(int j = 0; j < n; ++j)
        {
            const Scalar h_raw = step_scale * std::max(Scalar{1}, std::abs(x[j]));
            const Scalar temp = x[j] + h_raw;
            const Scalar h = temp - x[j];
            const Scalar orig = fd_x_buf_[j];

            fd_x_buf_[j] = orig + h;
            constraints(fd_x_buf_, c_plus);

            fd_x_buf_[j] = orig - h;
            constraints(fd_x_buf_, c_minus);

            J.col(j) = (c_plus - c_minus) / (Scalar{2} * h);

            fd_x_buf_[j] = orig;
        }
    }

    int n_eq{0};
    int n_ineq_upper{0};
    int n_ineq_lower{0};
    std::vector<int> eq_indices;
    std::vector<int> ineq_upper_indices;
    std::vector<int> ineq_lower_indices;

private:
    void eval_raw(const Eigen::VectorX<Scalar>& x) const
    {
        base::problem_->constraints(
            std::span<const Scalar>{x.data(), static_cast<std::size_t>(x.size())},
            std::span<Scalar>{raw_buf_.data(), raw_buf_.size()});
    }

    mutable std::vector<Scalar> raw_buf_;
    mutable Eigen::VectorX<Scalar> fd_x_buf_;
};

}

#endif
