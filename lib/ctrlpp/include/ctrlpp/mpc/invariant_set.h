#ifndef HPP_GUARD_CTRLPP_MPC_INVARIANT_SET_H
#define HPP_GUARD_CTRLPP_MPC_INVARIANT_SET_H

/// @brief Ellipsoidal and polytopic invariant set computation for terminal MPC constraints.
///
/// @cite mayne2000 -- Mayne et al., "Constrained model predictive control: Stability and optimality", 2000
/// @cite rawlings2017 -- Rawlings et al., "Model Predictive Control: Theory, Computation, and Design", 2017
/// @cite blanchini1999 -- Blanchini, "Set Invariance in Control", Automatica 35(11), 1999 (canonical survey on invariant sets)
/// @cite borrelli2017 -- Borrelli, Bemporad & Morari, "Predictive Control for Linear and Hybrid Systems", 2017, Ch. 10 (invariant-set computation, pre-image recursion)

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/control/dare.h"

#include "ctrlpp/mpc/terminal_set.h"

#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <vector>
#include <cstddef>
#include <numeric>
#include <optional>
#include <algorithm>

namespace ctrlpp
{

/// @cite mayne2000 -- Ellipsoidal invariant set from DARE solution, input, and state constraints
///
/// Caps alpha so the ellipsoid {x : x^T P x <= alpha} respects both the input limits (via the
/// LQR law u = K x) and the state box [x_min, x_max]. Requires u = 0 to be strictly interior to
/// [u_min, u_max] per component; returns terminal_set_error::input_zero_not_interior otherwise.
/// If the capped alpha is non-positive or non-finite (no face constrained it), the terminal set
/// is empty/unbounded and terminal_set_error::empty_terminal_set is returned -- alpha is never
/// silently shrunk and a non-finite alpha is never propagated.
template <typename Scalar, std::size_t NX, std::size_t NU>
auto compute_ellipsoidal_set(const Matrix<Scalar, NX, NX>& P,
                             const Matrix<Scalar, NU, NX>& K,
                             const Vector<Scalar, NU>& u_min,
                             const Vector<Scalar, NU>& u_max,
                             const Vector<Scalar, NX>& x_min,
                             const Vector<Scalar, NX>& x_max) -> ctrlpp::expected<ellipsoidal_set<Scalar, NX>, terminal_set_error>
{
    constexpr int nu = static_cast<int>(NU);
    constexpr int nx = static_cast<int>(NX);

    auto ldlt = P.ldlt();
    Scalar alpha = std::numeric_limits<Scalar>::infinity();

    // Input-face cap: alpha <= min(u_max(i)^2, u_min(i)^2) / (k_i^T P^{-1} k_i). The bound is
    // valid only when u = 0 lies strictly inside [u_min, u_max] for every component.
    for(int i = 0; i < nu; ++i)
    {
        if(!(u_min(i) < Scalar{0} && Scalar{0} < u_max(i)))
            return ctrlpp::unexpected(terminal_set_error::input_zero_not_interior);

        Vector<Scalar, NX> ki = K.row(i).transpose();
        auto Pinv_ki = ldlt.solve(ki).eval();

        Scalar denom = ki.dot(Pinv_ki);
        if(denom <= Scalar{0})
            continue;

        Scalar u_bound = std::min(u_max(i) * u_max(i), u_min(i) * u_min(i));
        alpha = std::min(alpha, u_bound / denom);
    }

    // State-face cap: the ellipsoid support along axis e_j is sqrt(alpha * e_j^T P^{-1} e_j);
    // keeping it inside [x_min(j), x_max(j)] requires alpha <= min(x_max(j)^2, x_min(j)^2) /
    // (e_j^T P^{-1} e_j). Infinite (unconstrained) state bounds leave alpha unchanged.
    for(int j = 0; j < nx; ++j)
    {
        Vector<Scalar, NX> ej = Vector<Scalar, NX>::Unit(j);
        auto Pinv_ej = ldlt.solve(ej).eval();

        Scalar denom = ej.dot(Pinv_ej);
        if(denom <= Scalar{0})
            continue;

        Scalar x_bound = std::min(x_max(j) * x_max(j), x_min(j) * x_min(j));
        alpha = std::min(alpha, x_bound / denom);
    }

    // alpha == 0 is the exact empty-set boundary; alpha == +inf means no face bounded the
    // ellipsoid at all. Either way there is no consistent, bounded terminal set,
    // so return an error rather than silently shrinking the set.
    if(!std::isfinite(alpha) || alpha <= Scalar{0})
        return ctrlpp::unexpected(terminal_set_error::empty_terminal_set);

    return ellipsoidal_set<Scalar, NX>{.P = P, .alpha = alpha};
}

namespace detail
{

/// Try to find a vertex at the intersection of the given hyperplane subset.
/// Returns std::nullopt if the system is singular or the vertex is infeasible.
template <typename Scalar, std::size_t NU>
auto try_vertex_at_intersection(const polytopic_set<Scalar, NU>& constraints,
                                const std::vector<int>& indices,
                                Scalar tol) -> std::optional<Vector<Scalar, NU>>
{
    constexpr int nu = static_cast<int>(NU);
    int m = static_cast<int>(constraints.H.rows());

    Eigen::Matrix<Scalar, nu, nu> H_sub;
    Eigen::Matrix<Scalar, nu, 1> h_sub;
    for(int i = 0; i < nu; ++i)
    {
        H_sub.row(i) = constraints.H.row(indices[static_cast<std::size_t>(i)]);
        h_sub(i) = constraints.h(indices[static_cast<std::size_t>(i)]);
    }

    auto qr = H_sub.colPivHouseholderQr();
    if(!qr.isInvertible())
        return std::nullopt;

    Vector<Scalar, NU> v = qr.solve(h_sub);

    for(int j = 0; j < m; ++j)
        if(constraints.H.row(j).dot(v) > constraints.h(j) + tol)
            return std::nullopt;

    return v;
}

/// Enumerate vertices of a polytopic constraint set via hyperplane intersection.
template <typename Scalar, std::size_t NU>
auto enumerate_polytope_vertices(const polytopic_set<Scalar, NU>& constraints, Scalar tol) -> std::vector<Vector<Scalar, NU>>
{
    constexpr int nu = static_cast<int>(NU);
    std::vector<Vector<Scalar, NU>> vertices;
    int m = static_cast<int>(constraints.H.rows());

    if(m < nu)
        return vertices;

    std::vector<int> indices(static_cast<std::size_t>(nu));
    std::iota(indices.begin(), indices.end(), 0);

    auto next_combination = [&]() -> bool
    {
        for(int i = nu - 1; i >= 0; --i)
        {
            auto idx = static_cast<std::size_t>(i);
            if(indices[idx] < m - nu + i)
            {
                ++indices[idx];
                for(int j = i + 1; j < nu; ++j)
                    indices[static_cast<std::size_t>(j)] = indices[static_cast<std::size_t>(j - 1)] + 1;
                return true;
            }
        }
        return false;
    };

    do
    {
        if(auto v = try_vertex_at_intersection<Scalar, NU>(constraints, indices, tol))
            vertices.push_back(*v);
    } while(next_combination());

    return vertices;
}

/// Compute the ROBUST (for-all-u) pre-image of a polytopic set under affine dynamics.
///
/// The pre-image constraints are stacked over ALL input vertices, so the resulting set is the
/// intersection { x : forall u in U, A x + B u in target } -- a point survives only if every
/// admissible input keeps it in the target. This yields a robust control-invariant set, not the
/// maximal one. The existential (exists-u) maximal-set variant is a separate operator and is not
/// implemented here.
template <typename Scalar, std::size_t NX, std::size_t NU>
auto compute_robust_pre_image(const Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>& H_curr,
                              const Eigen::VectorX<Scalar>& h_curr,
                              const Matrix<Scalar, NX, NX>& A_sys,
                              const Matrix<Scalar, NX, NU>& B_sys,
                              const std::vector<Vector<Scalar, NU>>& u_vertices)
    -> std::pair<Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>, Eigen::VectorX<Scalar>>
{
    constexpr int nx = static_cast<int>(NX);
    auto H_A = (H_curr * A_sys).eval();
    int n_faces = static_cast<int>(H_curr.rows());
    int n_verts = static_cast<int>(u_vertices.size());

    Eigen::Matrix<Scalar, Eigen::Dynamic, nx> H_pre(n_faces * n_verts, nx);
    Eigen::VectorX<Scalar> h_pre(n_faces * n_verts);

    for(int v = 0; v < n_verts; ++v)
    {
        auto Bu = (B_sys * u_vertices[static_cast<std::size_t>(v)]).eval();
        auto offset = (H_curr * Bu).eval();

        H_pre.middleRows(v * n_faces, n_faces) = H_A;
        h_pre.segment(v * n_faces, n_faces) = h_curr - offset;
    }

    return {H_pre, h_pre};
}

/// Check whether halfplane i is redundant w.r.t. already-kept halfplanes.
template <typename Scalar, std::size_t NX>
auto is_halfplane_redundant(const Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>& H,
                            const Eigen::VectorX<Scalar>& h,
                            int i,
                            const std::vector<int>& keep,
                            Scalar convergence_tol) -> bool
{
    Scalar row_norm = H.row(i).norm();
    Scalar h_normalized = h(i) / row_norm;
    Vector<Scalar, NX> dir_i = H.row(i).transpose() / row_norm;

    for(int j : keep)
    {
        Scalar row_norm_j = H.row(j).norm();
        Vector<Scalar, NX> dir_j = H.row(j).transpose() / row_norm_j;
        if(dir_i.dot(dir_j) > Scalar{1} - Scalar{1e-8})
        {
            Scalar h_norm_j = h(j) / row_norm_j;
            if(h_normalized >= h_norm_j - convergence_tol)
                return true;
        }
    }
    return false;
}

/// Extract a subset of rows from an H-representation.
template <typename Scalar, std::size_t NX>
auto extract_rows(const Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>& H,
                  const Eigen::VectorX<Scalar>& h,
                  const std::vector<int>& keep) -> std::pair<Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>, Eigen::VectorX<Scalar>>
{
    constexpr int nx = static_cast<int>(NX);
    Eigen::Matrix<Scalar, Eigen::Dynamic, nx> H_out(static_cast<int>(keep.size()), nx);
    Eigen::VectorX<Scalar> h_out(static_cast<int>(keep.size()));
    for(int i = 0; i < static_cast<int>(keep.size()); ++i)
    {
        H_out.row(i) = H.row(keep[static_cast<std::size_t>(i)]);
        h_out(i) = h(keep[static_cast<std::size_t>(i)]);
    }
    return {H_out, h_out};
}

/// Result of redundant-halfplane filtering: the reduced H-representation plus a flag reporting
/// whether the resource cap forced a truncation (so the reduced set is only an outer approximation).
template <typename Scalar, std::size_t NX>
struct filtered_halfplanes
{
    Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)> H;
    Eigen::VectorX<Scalar> h;
    bool truncated;
};

/// Remove redundant halfplanes from a polytopic H-representation.
///
/// The kept-halfplane count is bounded by a memory resource guard; if that bound is exceeded the
/// representation is truncated and `truncated` is set so callers can reject the (now unsound) set
/// instead of silently accepting a partial description.
template <typename Scalar, std::size_t NX>
auto filter_redundant_halfplanes(const Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>& H,
                                 const Eigen::VectorX<Scalar>& h,
                                 Scalar convergence_tol) -> filtered_halfplanes<Scalar, NX>
{
    int n_rows = static_cast<int>(H.rows());

    std::vector<int> keep;
    keep.reserve(static_cast<std::size_t>(n_rows));

    for(int i = 0; i < n_rows; ++i)
    {
        if(H.row(i).norm() < Scalar{1e-14})
            continue;
        if(!is_halfplane_redundant<Scalar, NX>(H, h, i, keep, convergence_tol))
            keep.push_back(i);
    }

    // Bounded-memory resource guard on the halfplane count. Hitting it truncates the H-representation
    // to an outer approximation, which is signaled via `truncated` rather than dropped silently.
    constexpr int max_halfplanes = 500;
    bool truncated = static_cast<int>(keep.size()) > max_halfplanes;
    if(truncated)
        keep.resize(static_cast<std::size_t>(max_halfplanes));

    auto [H_out, h_out] = extract_rows<Scalar, NX>(H, h, keep);
    return {.H = std::move(H_out), .h = std::move(h_out), .truncated = truncated};
}

}

/// Perform one backward-reachability iteration: compute the robust pre-image, merge, and filter.
/// Returns true if the set has converged (no new halfplanes added). Sets `truncated` if the filter
/// hit its resource cap during this step, so the caller can reject the resulting outer approximation.
template <typename Scalar, std::size_t NX, std::size_t NU>
auto backward_reachability_step(Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>& H_curr,
                                Eigen::VectorX<Scalar>& h_curr,
                                const Matrix<Scalar, NX, NX>& A_sys,
                                const Matrix<Scalar, NX, NU>& B_sys,
                                const std::vector<Vector<Scalar, NU>>& u_vertices,
                                Scalar convergence_tol,
                                bool& truncated) -> bool
{
    constexpr int nx = static_cast<int>(NX);
    auto [H_pre, h_pre] = detail::compute_robust_pre_image<Scalar, NX, NU>(H_curr, h_curr, A_sys, B_sys, u_vertices);

    int old_rows = static_cast<int>(H_curr.rows());
    int pre_rows = static_cast<int>(H_pre.rows());

    Eigen::Matrix<Scalar, Eigen::Dynamic, nx> H_next(old_rows + pre_rows, nx);
    Eigen::VectorX<Scalar> h_next(old_rows + pre_rows);
    H_next.topRows(old_rows) = H_curr;
    h_next.head(old_rows) = h_curr;
    H_next.bottomRows(pre_rows) = H_pre;
    h_next.tail(pre_rows) = h_pre;

    auto filtered = detail::filter_redundant_halfplanes<Scalar, NX>(H_next, h_next, convergence_tol);
    truncated = filtered.truncated;
    bool converged = (static_cast<int>(filtered.H.rows()) == old_rows);

    H_curr = std::move(filtered.H);
    h_curr = std::move(filtered.h);
    return converged;
}

/// Compute the polytopic ROBUST control-invariant set via backward reachability.
///
/// Each iteration intersects the robust (for-all-u) pre-image with the current set, so the
/// fixed point is the robust control-invariant set (the states kept inside the state box under
/// every admissible input), not the maximal control-invariant set. Returns
/// terminal_set_error::halfplanes_truncated if the halfplane resource cap fires (the set is only
/// an outer approximation), terminal_set_error::not_converged if the iteration budget is exhausted
/// before a fixed point, and terminal_set_error::empty_terminal_set if the input polytope has no
/// enumerable vertices. On success the returned set is a converged robust control-invariant set.
/// @cite rawlings2017 -- Ch. 2 (invariant set computation)
template <typename Scalar, std::size_t NX, std::size_t NU>
auto compute_polytopic_invariant_set(const Matrix<Scalar, NX, NX>& A_sys,
                                     const Matrix<Scalar, NX, NU>& B_sys,
                                     const polytopic_set<Scalar, NX>& state_constraints,
                                     const polytopic_set<Scalar, NU>& input_constraints,
                                     int max_iterations = 100,
                                     Scalar convergence_tol = Scalar{1e-6}) -> ctrlpp::expected<polytopic_set<Scalar, NX>, terminal_set_error>
{
    static_assert(NX <= 4, "Polytopic invariant set computation restricted to NX <= 4");

    auto u_vertices = detail::enumerate_polytope_vertices<Scalar, NU>(input_constraints, convergence_tol);
    if(u_vertices.empty())
        return ctrlpp::unexpected(terminal_set_error::empty_terminal_set);

    auto H_curr = state_constraints.H;
    auto h_curr = state_constraints.h;

    for(int iter = 0; iter < max_iterations; ++iter)
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
};

/// @cite mayne2000 -- Terminal cost and constraint set from DARE + LQR
///
/// Returns terminal_set_error::dare_failed if the Riccati solve fails, or propagates the
/// terminal_set_error from compute_ellipsoidal_set (input-sign / empty-set) on the error channel.
template <typename Scalar, std::size_t NX, std::size_t NU>
auto terminal_ingredients(
    const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NX, NU>& B, const Matrix<Scalar, NX, NX>& Q, const Matrix<Scalar, NU, NU>& R,
    const Vector<Scalar, NU>& u_min, const Vector<Scalar, NU>& u_max, const Vector<Scalar, NX>& x_min, const Vector<Scalar, NX>& x_max)
    -> ctrlpp::expected<terminal_ingredients_result<Scalar, NX, NU>, terminal_set_error>
{
    auto P_result = dare<Scalar, NX, NU>(A, B, Q, R);
    if(!P_result)
        return ctrlpp::unexpected(terminal_set_error::dare_failed);

    auto P = P_result->P;
    Matrix<Scalar, NU, NU> RpBtPB = R + B.transpose() * P * B;
    Matrix<Scalar, NU, NX> K = -(RpBtPB.ldlt().solve(B.transpose() * P * A));
    auto eset = compute_ellipsoidal_set<Scalar, NX, NU>(P, K, u_min, u_max, x_min, x_max);
    if(!eset)
        return ctrlpp::unexpected(eset.error());

    return terminal_ingredients_result<Scalar, NX, NU>{.Qf = P, .set = *eset};
}

}

#endif
