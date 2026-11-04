#ifndef HPP_GUARD_CTRLPP_MPC_INVARIANT_SET_H
#define HPP_GUARD_CTRLPP_MPC_INVARIANT_SET_H

// @cite mayne2000 -- Mayne et al., "Constrained model predictive control: Stability and optimality", 2000
// @cite rawlings2017 -- Rawlings et al., "Model Predictive Control: Theory, Computation, and Design", 2017
// @cite blanchini1999 -- Blanchini, "Set Invariance in Control", Automatica 35(11), 1999 (canonical survey on invariant sets)
// @cite borrelli2017 -- Borrelli, Bemporad & Morari, "Predictive Control for Linear and Hybrid Systems", 2017, Ch. 10 (invariant-set computation, pre-image recursion)

#include "ctrlpp/mpc/terminal_set.h"
#include "ctrlpp/mpc/terminal_cost.h"
#include "ctrlpp/mpc/output_admissible_set.h"

#include "ctrlpp/mpc/detail/polytope_vertices.h"
#include "ctrlpp/mpc/detail/invariant_set_rows.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/control/dare.h"

#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <algorithm>

namespace ctrlpp
{

// @cite mayne2000 -- Ellipsoidal invariant set from DARE solution, input, and state constraints
// The ellipsoid's support along d is sqrt(alpha d^T P^{-1} d), so a face [lower, upper] along d caps alpha at
// min(upper^2, lower^2) / (d^T P^{-1} d), with d = k_i for an input row and d = e_j for a state axis. The input caps
// hold only when u = 0 lies strictly inside [u_min, u_max]; an infinite state bound leaves alpha unchanged. alpha == 0
// is the exact empty-set boundary, and alpha == +inf means no face bounded the ellipsoid.
template <typename Scalar, std::size_t NX, std::size_t NU>
ctrlpp::expected<ellipsoidal_set<Scalar, NX>, terminal_set_error> compute_ellipsoidal_set(const Matrix<Scalar, NX, NX>& P, const Matrix<Scalar, NU, NX>& K,
        const Vector<Scalar, NU>& u_min, const Vector<Scalar, NU>& u_max, const Vector<Scalar, NX>& x_min, const Vector<Scalar, NX>& x_max)
{
    auto ldlt = P.ldlt();
    Scalar alpha = std::numeric_limits<Scalar>::infinity();
    for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NU); ++i)
    {
        if(!(u_min(i) < Scalar{0} && Scalar{0} < u_max(i)))
            return ctrlpp::unexpected(terminal_set_error::input_zero_not_interior);
        Vector<Scalar, NX> ki = K.row(i).transpose();
        Scalar denom = ki.dot(ldlt.solve(ki).eval());
        if(denom <= Scalar{0})
            continue;
        alpha = std::min(alpha, std::min(u_max(i) * u_max(i), u_min(i) * u_min(i)) / denom);
    }
    for(Eigen::Index j = 0; j < static_cast<Eigen::Index>(NX); ++j)
    {
        Vector<Scalar, NX> ej = Vector<Scalar, NX>::Unit(j);
        auto Pinv_ej = ldlt.solve(ej).eval();
        Scalar denom = ej.dot(Pinv_ej);
        if(denom <= Scalar{0})
            continue;
        alpha = std::min(alpha, std::min(x_max(j) * x_max(j), x_min(j) * x_min(j)) / denom);
    }
    if(!std::isfinite(alpha) || alpha <= Scalar{0})
        return ctrlpp::unexpected(terminal_set_error::empty_terminal_set);
    return ellipsoidal_set<Scalar, NX>{.P = P, .alpha = alpha};
}

// Converged when filtering leaves the row count unchanged; `truncated` reports that the filter hit its cap.
template <typename Scalar, std::size_t NX, std::size_t NU>
bool backward_reachability_step(Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>& H_curr, Eigen::VectorX<Scalar>& h_curr, const Matrix<Scalar, NX, NX>& A_sys,
        const Matrix<Scalar, NX, NU>& B_sys, const std::vector<Vector<Scalar, NU>>& u_vertices, Scalar convergence_tol, bool& truncated)
{
    constexpr int nx = static_cast<int>(NX);
    auto [H_pre, h_pre] = detail::compute_robust_pre_image<Scalar, NX, NU>(H_curr, h_curr, A_sys, B_sys, u_vertices);

    Eigen::Index old_rows = H_curr.rows();
    Eigen::Index pre_rows = H_pre.rows();

    Eigen::Matrix<Scalar, Eigen::Dynamic, nx> H_next(old_rows + pre_rows, nx);
    Eigen::VectorX<Scalar> h_next(old_rows + pre_rows);
    H_next.topRows(old_rows) = H_curr;
    h_next.head(old_rows) = h_curr;
    H_next.bottomRows(pre_rows) = H_pre;
    h_next.tail(pre_rows) = h_pre;

    auto filtered = detail::filter_redundant_halfplanes<Scalar, NX>(H_next, h_next, convergence_tol);
    truncated = filtered.truncated;
    bool converged = (filtered.H.rows() == old_rows);

    H_curr = std::move(filtered.H);
    h_curr = std::move(filtered.h);
    return converged;
}

// @cite rawlings2017 -- Ch. 2 (invariant set computation)
template <typename Scalar, std::size_t NX, std::size_t NU>
ctrlpp::expected<polytopic_set<Scalar, NX>, terminal_set_error> compute_polytopic_invariant_set(const Matrix<Scalar, NX, NX>& A_sys, const Matrix<Scalar, NX, NU>& B_sys,
        const polytopic_set<Scalar, NX>& state_constraints, const polytopic_set<Scalar, NU>& input_constraints, std::int32_t max_iterations = 100, Scalar convergence_tol = Scalar{1e-6})
{
    static_assert(NX <= 4, "Polytopic invariant set computation restricted to NX <= 4");

    auto u_vertices = detail::enumerate_polytope_vertices<Scalar, NU>(input_constraints, convergence_tol);
    if(u_vertices.empty())
        return ctrlpp::unexpected(terminal_set_error::empty_terminal_set);

    auto H_curr = state_constraints.H;
    auto h_curr = state_constraints.h;

    for(std::int32_t iter = 0; iter < max_iterations; ++iter)
    {
        bool truncated = false;
        bool converged = backward_reachability_step<Scalar, NX, NU>(H_curr, h_curr, A_sys, B_sys, u_vertices, convergence_tol, truncated);
        if(truncated)
            return ctrlpp::unexpected(terminal_set_error::halfplanes_truncated);
        if(converged)
            return polytopic_set<Scalar, NX>{.H = std::move(H_curr), .h = std::move(h_curr)};
    }

    return ctrlpp::unexpected(terminal_set_error::not_converged);
}

template <typename Scalar, std::size_t NX, std::size_t NU>
struct terminal_ingredients_result
{
    Matrix<Scalar, NX, NX> Qf;
    ellipsoidal_set<Scalar, NX> set;
    ctrlpp::expected<polytopic_set<Scalar, NX>, terminal_set_error> invariant_set;
};

// @cite mayne2000 -- Terminal cost and constraint sets from DARE + LQR
template <typename Scalar, std::size_t NX, std::size_t NU>
ctrlpp::expected<terminal_ingredients_result<Scalar, NX, NU>, terminal_ingredients_error> terminal_ingredients(const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NX, NU>& B,
        const Matrix<Scalar, NX, NX>& Q, const Matrix<Scalar, NU, NU>& R, const Vector<Scalar, NU>& u_min, const Vector<Scalar, NU>& u_max, const Vector<Scalar, NX>& x_min,
        const Vector<Scalar, NX>& x_max)
{
    if(detail::riccati_input_weight_singular<Scalar, NX, NU>(B, Q, R))
        return ctrlpp::unexpected(terminal_ingredients_error{terminal_set_error::singular_input_weight});
    auto P_result = dare<Scalar, NX, NU>(A, B, Q, R);
    if(!P_result)
        return ctrlpp::unexpected(terminal_ingredients_error{P_result.error()});

    auto P = P_result->P;
    Matrix<Scalar, NU, NU> RpBtPB = R + B.transpose() * P * B;
    Matrix<Scalar, NU, NX> K = -(RpBtPB.ldlt().solve(B.transpose() * P * A));
    auto eset = compute_ellipsoidal_set<Scalar, NX, NU>(P, K, u_min, u_max, x_min, x_max);
    if(!eset)
        return ctrlpp::unexpected(terminal_ingredients_error{eset.error()});

    terminal_ingredients_result<Scalar, NX, NU> result{.Qf = P, .set = *eset, .invariant_set = ctrlpp::unexpected(terminal_set_error::dimension_unsupported)};
    if constexpr(NX <= 4)
        result.invariant_set = compute_output_admissible_set<Scalar, NX, NU>(A + B * K, K, u_min, u_max, x_min, x_max);
    return result;
}

}

#endif
