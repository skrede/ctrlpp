#ifndef HPP_GUARD_CTRLPP_MPC_INVARIANT_SET_H
#define HPP_GUARD_CTRLPP_MPC_INVARIANT_SET_H

/// @brief Ellipsoidal and polytopic invariant set computation for terminal MPC constraints.
///
/// @cite mayne2000 -- Mayne et al., "Constrained model predictive control: Stability and optimality", 2000
/// @cite rawlings2017 -- Rawlings et al., "Model Predictive Control: Theory, Computation, and Design", 2017

#include "ctrlpp/mpc/terminal_set.h"
#include "ctrlpp/control/dare.h"
#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <optional>
#include <vector>

namespace ctrlpp
{

/// @cite mayne2000 -- Ellipsoidal invariant set from DARE solution and input constraints
template <typename Scalar, std::size_t NX, std::size_t NU>
auto compute_ellipsoidal_set(const Matrix<Scalar, NX, NX>& P, const Matrix<Scalar, NU, NX>& K, const Vector<Scalar, NU>& u_min, const Vector<Scalar, NU>& u_max) -> ellipsoidal_set<Scalar, NX>
{
    constexpr int nu = static_cast<int>(NU);

    auto ldlt = P.ldlt();
    Scalar alpha = std::numeric_limits<Scalar>::infinity();

    for(int i = 0; i < nu; ++i)
    {
        Vector<Scalar, NX> ki = K.row(i).transpose();
        auto Pinv_ki = ldlt.solve(ki).eval();

        Scalar denom = ki.dot(Pinv_ki);
        if(denom <= Scalar{0})
            continue;

        Scalar u_bound = std::min(u_max(i) * u_max(i), u_min(i) * u_min(i));
        alpha = std::min(alpha, u_bound / denom);
    }

    return {.P = P, .alpha = alpha};
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

/// Compute pre-image of a polytopic set under affine dynamics.
template <typename Scalar, std::size_t NX, std::size_t NU>
auto compute_pre_image(const Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>& H_curr,
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

/// Remove redundant halfplanes from a polytopic H-representation.
template <typename Scalar, std::size_t NX>
auto filter_redundant_halfplanes(const Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>& H,
                                 const Eigen::VectorX<Scalar>& h,
                                 Scalar convergence_tol) -> std::pair<Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>, Eigen::VectorX<Scalar>>
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

    constexpr int max_halfplanes = 500;
    if(static_cast<int>(keep.size()) > max_halfplanes)
        keep.resize(static_cast<std::size_t>(max_halfplanes));

    return extract_rows<Scalar, NX>(H, h, keep);
}

}

/// Perform one backward-reachability iteration: compute pre-image, merge, and filter.
/// Returns true if the set has converged (no new halfplanes added).
template <typename Scalar, std::size_t NX, std::size_t NU>
auto backward_reachability_step(Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NX)>& H_curr,
                                Eigen::VectorX<Scalar>& h_curr,
                                const Matrix<Scalar, NX, NX>& A_sys,
                                const Matrix<Scalar, NX, NU>& B_sys,
                                const std::vector<Vector<Scalar, NU>>& u_vertices,
                                Scalar convergence_tol) -> bool
{
    constexpr int nx = static_cast<int>(NX);
    auto [H_pre, h_pre] = detail::compute_pre_image<Scalar, NX, NU>(H_curr, h_curr, A_sys, B_sys, u_vertices);

    int old_rows = static_cast<int>(H_curr.rows());
    int pre_rows = static_cast<int>(H_pre.rows());

    Eigen::Matrix<Scalar, Eigen::Dynamic, nx> H_next(old_rows + pre_rows, nx);
    Eigen::VectorX<Scalar> h_next(old_rows + pre_rows);
    H_next.topRows(old_rows) = H_curr;
    h_next.head(old_rows) = h_curr;
    H_next.bottomRows(pre_rows) = H_pre;
    h_next.tail(pre_rows) = h_pre;

    auto [H_filtered, h_filtered] = detail::filter_redundant_halfplanes<Scalar, NX>(H_next, h_next, convergence_tol);
    bool converged = (static_cast<int>(H_filtered.rows()) == old_rows);

    H_curr = std::move(H_filtered);
    h_curr = std::move(h_filtered);
    return converged;
}

/// Compute polytopic maximal control-invariant set via backward reachability.
/// @cite rawlings2017 -- Ch. 2 (invariant set computation)
template <typename Scalar, std::size_t NX, std::size_t NU>
auto compute_polytopic_invariant_set(const Matrix<Scalar, NX, NX>& A_sys,
                                     const Matrix<Scalar, NX, NU>& B_sys,
                                     const polytopic_set<Scalar, NX>& state_constraints,
                                     const polytopic_set<Scalar, NU>& input_constraints,
                                     int max_iterations = 100,
                                     Scalar convergence_tol = Scalar{1e-6}) -> std::optional<polytopic_set<Scalar, NX>>
{
    static_assert(NX <= 4, "Polytopic invariant set computation restricted to NX <= 4");

    auto u_vertices = detail::enumerate_polytope_vertices<Scalar, NU>(input_constraints, convergence_tol);
    if(u_vertices.empty())
        return std::nullopt;

    auto H_curr = state_constraints.H;
    auto h_curr = state_constraints.h;

    for(int iter = 0; iter < max_iterations; ++iter)
        if(backward_reachability_step<Scalar, NX, NU>(H_curr, h_curr, A_sys, B_sys, u_vertices, convergence_tol))
            return polytopic_set<Scalar, NX>{.H = std::move(H_curr), .h = std::move(h_curr)};

    return polytopic_set<Scalar, NX>{.H = std::move(H_curr), .h = std::move(h_curr)};
}

template <typename Scalar, std::size_t NX, std::size_t NU>
struct terminal_ingredients_result
{
    Matrix<Scalar, NX, NX> Qf;
    ellipsoidal_set<Scalar, NX> set;
};

/// @cite mayne2000 -- Terminal cost and constraint set from DARE + LQR
template <typename Scalar, std::size_t NX, std::size_t NU>
auto terminal_ingredients(
    const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NX, NU>& B, const Matrix<Scalar, NX, NX>& Q, const Matrix<Scalar, NU, NU>& R, const Vector<Scalar, NU>& u_min, const Vector<Scalar, NU>& u_max)
    -> std::optional<terminal_ingredients_result<Scalar, NX, NU>>
{
    auto P_opt = dare<Scalar, NX, NU>(A, B, Q, R);
    if(!P_opt)
        return std::nullopt;

    auto P = P_opt.value();
    Matrix<Scalar, NU, NU> RpBtPB = R + B.transpose() * P * B;
    Matrix<Scalar, NU, NX> K = -(RpBtPB.ldlt().solve(B.transpose() * P * A));
    auto eset = compute_ellipsoidal_set<Scalar, NX, NU>(P, K, u_min, u_max);

    return terminal_ingredients_result<Scalar, NX, NU>{.Qf = P, .set = eset};
}

}

#endif
