#ifndef HPP_GUARD_CTRLPP_MPC_ARGMIN_PROBLEM_H
#define HPP_GUARD_CTRLPP_MPC_ARGMIN_PROBLEM_H

#include "ctrlpp/mpc/nlp_solver.h"

#include <Eigen/Core>

#include <span>
#include <cmath>
#include <vector>
#include <cassert>
#include <cstddef>
#include <algorithm>
#include <type_traits>

namespace ctrlpp
{

namespace detail
{

// Optional compile-time constraint-count channel (argmin SEED-044). argmin
// detects a problem's `static constexpr int constraint_count` by PRESENCE
// (has_constraint_count<P>, argmin/formulation/concepts.h). Presence opts into a
// compile-time M for the QP result-multiplier storage; absence reproduces the
// dynamic-M floor. A -1 sentinel would NOT work as the dynamic signal: argmin
// reads any present value as the cap and asserts runtime_m <= cap, so an absent
// member -- not a sentinel value -- is how we request dynamic behavior.
template <int MaxM>
struct argmin_constraint_bound
{
    static constexpr int constraint_count = MaxM;
};

template <>
struct argmin_constraint_bound<Eigen::Dynamic>
{
    // Intentionally empty: no constraint_count member => argmin dynamic-M path.
};

}

// argmin_problem / argmin_constrained_problem bridge the ctrlpp NLP contract to
// argmin's problem interface. They are parameterized on two independent
// compile-time axes: the decision dimension NV and the constraint-count cap MaxM.
//
//   * NV == Eigen::Dynamic (the DEFAULT) reproduces the original runtime-erased
//     bridge byte-for-byte: it binds an nlp_problem<Scalar>, its decision-vector
//     storage is Eigen::VectorX<Scalar>, and `allocation_free` is false. Every
//     existing caller compiles unchanged.
//   * NV != Eigen::Dynamic sizes the decision-vector storage with fixed-size
//     Eigen types (Eigen::Vector<Scalar, NV>), binds an nlp_problem_static<Scalar,
//     NV>, and exposes `allocation_free == true`. This pins the DECISION axis.
//
// The constraint axis is bound separately via MaxM (argmin SEED-044):
//
//   * MaxM == Eigen::Dynamic (the DEFAULT) carries NO `constraint_count` member,
//     so argmin's has_constraint_count<P> is false and the constraint-axis result
//     multipliers stay dynamic-but-preallocated -- argmin's fixed-N NW-SQP floor,
//     where argmin heap-allocates the per-call multiplier storage.
//   * MaxM != Eigen::Dynamic exposes `static constexpr int constraint_count = MaxM`
//     (an UPPER BOUND; runtime_m <= MaxM), which argmin threads
//     constraint_count -> state_type<P>::M -> active_set_qp_solver<Scalar, NV, M>,
//     flipping bounds_fit_inline true so the result-multiplier storage is inline.
//     Binding BOTH axes (NV and MaxM) yields the strict-zero steady-state solve
//     (`strict_allocation_free`). Box bounds are free via argmin's +2N slack, so
//     MaxM counts equality + general-inequality rows only.
template <typename Scalar, int NV = Eigen::Dynamic, int MaxM = Eigen::Dynamic>
class argmin_problem : public detail::argmin_constraint_bound<MaxM>
{
public:
    static constexpr int problem_dimension = NV;

    // Compile-time marker: the fixed-NV specialization backs its decision-vector
    // storage with fixed-size Eigen types and does not allocate on the
    // solve/step path for the DECISION axis, whereas the dynamic default does.
    // Consumed by the static-path test's static_assert.
    static constexpr bool allocation_free = (NV != Eigen::Dynamic);

    // Compile-time marker: the constraint axis is bound to a compile-time cap, so
    // argmin's per-call result-multiplier storage is inline rather than heap.
    static constexpr bool constraint_bounded = (MaxM != Eigen::Dynamic);

    // Joint strict-zero marker: BOTH the decision axis (fixed-size storage) and
    // the constraint axis (inline multipliers) are bound, so the steady-state
    // solve/step path performs zero heap allocation end-to-end. Consumed by
    // nmpc_static_nomalloc_test's static_assert.
    static constexpr bool strict_allocation_free = allocation_free && constraint_bounded;

    // The bound problem type follows NV: the non-erased static contract on the
    // fixed-NV path, the runtime-erased contract on the dynamic default.
    using problem_type = std::conditional_t<NV == Eigen::Dynamic,
        nlp_problem<Scalar>,
        nlp_problem_static<Scalar, NV>>;

    // Decision-vector storage type: fixed-size for NV != Eigen::Dynamic,
    // identical to Eigen::VectorX<Scalar> when NV == Eigen::Dynamic.
    using decision_vector = Eigen::Vector<Scalar, NV>;

    void bind(const problem_type& prob)
    {
        problem_ = &prob;
    }

    auto value(const decision_vector& x) const -> Scalar
    {
        return problem_->cost(
            std::span<const Scalar>{x.data(), static_cast<std::size_t>(x.size())});
    }

    auto dimension() const -> int
    {
        return problem_->n_vars;
    }

    void gradient(const decision_vector& x, decision_vector& g) const
    {
        // Defensive size check: argmin owns the output buffer and must size it
        // to the problem dimension before this write. A debug-only assert
        // (compiled out under NDEBUG, no allocation, no throw) guards the hot
        // path without disturbing the RT/no-exceptions posture.
        assert(g.size() == static_cast<Eigen::Index>(dimension()));

        problem_->gradient(
            std::span<const Scalar>{x.data(), static_cast<std::size_t>(x.size())},
            std::span<Scalar>{g.data(), static_cast<std::size_t>(g.size())});
    }

    auto lower_bounds() const -> decision_vector
    {
        return problem_->x_lower;
    }

    auto upper_bounds() const -> decision_vector
    {
        return problem_->x_upper;
    }

protected:
    const problem_type* problem_{nullptr};
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

    void partition(const problem_type& prob)
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

        // All per-call workspaces are sized ONCE here (at setup), so the
        // solve/step path that follows performs no heap allocation. J_raw_ is
        // the analytic-Jacobian reorder buffer hoisted out of constraint_jacobian
        // (it was a per-call temporary; on the static path a per-call heap
        // allocation would break the allocation-free contract).
        raw_buf_.resize(static_cast<std::size_t>(prob.n_constraints));
        fd_x_buf_.resize(prob.n_vars);
        c_plus_.resize(n_eq + n_ineq_upper + n_ineq_lower);
        c_minus_.resize(n_eq + n_ineq_upper + n_ineq_lower);
        J_raw_.resize(prob.n_constraints, prob.n_vars);
    }

    void constraints(const decision_vector& x, Eigen::VectorX<Scalar>& c_out) const
    {
        // Defensive size check: c_out is caller/argmin-provided and must hold
        // exactly the partitioned constraint count (equalities + upper/lower
        // inequalities) before any write below. A debug-only assert (compiled
        // out under NDEBUG, no allocation, no throw) prevents an out-of-bounds
        // write on the hot path without disturbing the RT/no-exceptions posture.
        assert(c_out.size()
            == static_cast<Eigen::Index>(n_eq + n_ineq_upper + n_ineq_lower));

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

    void constraint_jacobian(const decision_vector& x, Eigen::MatrixX<Scalar>& J) const
    {
        const int m = n_eq + n_ineq_upper + n_ineq_lower;
        const int n = base::problem_->n_vars;
        J.resize(m, n);

        if(base::problem_->constraint_jacobian)
        {
            // J_raw_ is a pre-sized member (partition() sized it to
            // n_constraints x n_vars), so this reorder path allocates nothing.
            base::problem_->constraint_jacobian(
                std::span<const Scalar>{x.data(), static_cast<std::size_t>(x.size())},
                std::span<Scalar>{J_raw_.data(), static_cast<std::size_t>(J_raw_.size())});

            int idx = 0;
            for(auto i : eq_indices)
                J.row(idx++) = J_raw_.row(i);
            for(auto i : ineq_upper_indices)
                J.row(idx++) = J_raw_.row(i);
            for(auto i : ineq_lower_indices)
                J.row(idx++) = -J_raw_.row(i);

            return;
        }

        const auto step_scale = std::cbrt(std::numeric_limits<Scalar>::epsilon());

        fd_x_buf_ = x;

        for(int j = 0; j < n; ++j)
        {
            const Scalar h_raw = step_scale * std::max(Scalar{1}, std::abs(x[j]));
            const Scalar temp = x[j] + h_raw;
            const Scalar h = temp - x[j];
            const Scalar orig = fd_x_buf_[j];

            fd_x_buf_[j] = orig + h;
            constraints(fd_x_buf_, c_plus_);

            fd_x_buf_[j] = orig - h;
            constraints(fd_x_buf_, c_minus_);

            J.col(j) = (c_plus_ - c_minus_) / (Scalar{2} * h);

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
    void eval_raw(const decision_vector& x) const
    {
        base::problem_->constraints(
            std::span<const Scalar>{x.data(), static_cast<std::size_t>(x.size())},
            std::span<Scalar>{raw_buf_.data(), raw_buf_.size()});
    }

    mutable std::vector<Scalar> raw_buf_;
    mutable decision_vector fd_x_buf_;
    mutable Eigen::VectorX<Scalar> c_plus_;
    mutable Eigen::VectorX<Scalar> c_minus_;
    mutable Eigen::MatrixX<Scalar> J_raw_;
};

}

#endif
