#ifndef HPP_GUARD_CTRLPP_MPC_ARGMIN_PROBLEM_H
#define HPP_GUARD_CTRLPP_MPC_ARGMIN_PROBLEM_H

#include "ctrlpp/mpc/nlp_solver.h"

#include "ctrlpp/mpc/detail/argmin_problem_rows.h"

#include <Eigen/Core>

#include <span>
#include <cmath>
#include <limits>
#include <vector>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <type_traits>

namespace ctrlpp
{

template <typename Scalar, int NV = Eigen::Dynamic, int MaxM = Eigen::Dynamic>
class argmin_problem : public detail::argmin_constraint_bound<MaxM>
{
public:
    static constexpr int problem_dimension = NV;
    static constexpr bool allocation_free = (NV != Eigen::Dynamic);
    static constexpr bool constraint_bounded = (MaxM != Eigen::Dynamic);
    static constexpr bool strict_allocation_free = allocation_free && constraint_bounded;

    using problem_type = std::conditional_t<NV == Eigen::Dynamic, nlp_problem<Scalar>, nlp_problem_static<Scalar, NV>>;
    using decision_vector = Eigen::Vector<Scalar, NV>;

    argmin_problem()
        : m_problem{nullptr}
    {
    }

    void bind(const problem_type& prob) { m_problem = &prob; }

    Scalar value(const decision_vector& x) const { return problem().cost_cb(std::span<const Scalar>{x.data(), static_cast<std::size_t>(x.size())}); }

    std::int32_t dimension() const { return static_cast<std::int32_t>(problem().n_vars); }

    void gradient(const decision_vector& x, decision_vector& g) const
    {
        // argmin owns g and sizes it to the problem dimension before this call.
        assert(g.size() == static_cast<Eigen::Index>(dimension()));

        problem().gradient_cb(std::span<const Scalar>{x.data(), static_cast<std::size_t>(x.size())}, std::span<Scalar>{g.data(), static_cast<std::size_t>(g.size())});
    }

    decision_vector lower_bounds() const { return problem().x_lower; }

    decision_vector upper_bounds() const { return problem().x_upper; }

protected:
    // Null until bind(). An optional reference here keeps GCC from inlining value() and gradient() into argmin's step, which grows its frame.
    const problem_type* m_problem;

    const problem_type& problem() const { return *m_problem; }
};

template <typename Scalar, int NV = Eigen::Dynamic, int MaxM = Eigen::Dynamic>
class argmin_constrained_problem : public argmin_problem<Scalar, NV, MaxM>
{
    using base = argmin_problem<Scalar, NV, MaxM>;

public:
    using base::problem_dimension;
    using base::allocation_free;
    using base::constraint_bounded;
    using base::strict_allocation_free;
    using problem_type = typename base::problem_type;
    using decision_vector = typename base::decision_vector;

    argmin_constrained_problem()
        : n_eq{0}
        , n_ineq_upper{0}
        , n_ineq_lower{0}
    {
    }

    void partition(const problem_type& prob)
    {
        base::bind(prob);

        detail::classify_constraint_rows(prob, eq_indices, ineq_upper_indices, ineq_lower_indices);

        n_eq = static_cast<std::int32_t>(eq_indices.size());
        n_ineq_upper = static_cast<std::int32_t>(ineq_upper_indices.size());
        n_ineq_lower = static_cast<std::int32_t>(ineq_lower_indices.size());

        // Every per-call workspace is sized here, so the solve and step path allocates nothing.
        m_raw_buf.resize(static_cast<std::size_t>(prob.n_constraints));
        m_fd_x_buf.resize(prob.n_vars);
        m_c_plus.resize(n_eq + n_ineq_upper + n_ineq_lower);
        m_c_minus.resize(n_eq + n_ineq_upper + n_ineq_lower);
        m_J_raw.resize(prob.n_constraints, prob.n_vars);
    }

    void constraints(const decision_vector& x, Eigen::VectorX<Scalar>& c_out) const
    {
        assert(c_out.size() == static_cast<Eigen::Index>(n_eq + n_ineq_upper + n_ineq_lower));

        eval_raw(x);

        // argmin takes an inequality row as satisfied when it is non-negative
        // (argmin/test_functions/hock_schittkowski.h, "Inequality convention").
        Eigen::Index idx = 0;

        for(auto i : eq_indices)
            c_out[idx++] = m_raw_buf[static_cast<std::size_t>(i)] - base::problem().c_lower[i];

        for(auto i : ineq_upper_indices)
            c_out[idx++] = base::problem().c_upper[i] - m_raw_buf[static_cast<std::size_t>(i)];

        for(auto i : ineq_lower_indices)
            c_out[idx++] = m_raw_buf[static_cast<std::size_t>(i)] - base::problem().c_lower[i];
    }

    std::int32_t num_equality() const { return n_eq; }

    std::int32_t num_inequality() const { return n_ineq_upper + n_ineq_lower; }

    void constraint_jacobian(const decision_vector& x, Eigen::MatrixX<Scalar>& J) const
    {
        const Eigen::Index m = n_eq + n_ineq_upper + n_ineq_lower;
        const Eigen::Index n = base::problem().n_vars;
        J.resize(m, n);

        if(base::problem().constraint_jacobian_cb)
        {
            base::problem().constraint_jacobian_cb(
                    std::span<const Scalar>{x.data(), static_cast<std::size_t>(x.size())}, std::span<Scalar>{m_J_raw.data(), static_cast<std::size_t>(m_J_raw.size())});

            Eigen::Index idx = 0;
            for(auto i : eq_indices)
                J.row(idx++) = m_J_raw.row(i);
            for(auto i : ineq_upper_indices)
                J.row(idx++) = -m_J_raw.row(i);
            for(auto i : ineq_lower_indices)
                J.row(idx++) = m_J_raw.row(i);

            return;
        }

        const auto step_scale = std::cbrt(std::numeric_limits<Scalar>::epsilon());

        m_fd_x_buf = x;

        for(Eigen::Index j = 0; j < n; ++j)
        {
            const Scalar h_raw = step_scale * std::max(Scalar{1}, std::abs(x[j]));
            const Scalar temp = x[j] + h_raw;
            const Scalar h = temp - x[j];
            const Scalar orig = m_fd_x_buf[j];

            m_fd_x_buf[j] = orig + h;
            constraints(m_fd_x_buf, m_c_plus);

            m_fd_x_buf[j] = orig - h;
            constraints(m_fd_x_buf, m_c_minus);

            J.col(j) = (m_c_plus - m_c_minus) / (Scalar{2} * h);

            m_fd_x_buf[j] = orig;
        }
    }

    std::int32_t n_eq;
    std::int32_t n_ineq_upper;
    std::int32_t n_ineq_lower;
    std::vector<Eigen::Index> eq_indices;
    std::vector<Eigen::Index> ineq_upper_indices;
    std::vector<Eigen::Index> ineq_lower_indices;

private:
    mutable decision_vector m_fd_x_buf;
    mutable std::vector<Scalar> m_raw_buf;
    mutable Eigen::MatrixX<Scalar> m_J_raw;
    mutable Eigen::VectorX<Scalar> m_c_plus;
    mutable Eigen::VectorX<Scalar> m_c_minus;

    void eval_raw(const decision_vector& x) const
    {
        base::problem().constraints_cb(std::span<const Scalar>{x.data(), static_cast<std::size_t>(x.size())}, std::span<Scalar>{m_raw_buf.data(), m_raw_buf.size()});
    }
};

}

#endif
