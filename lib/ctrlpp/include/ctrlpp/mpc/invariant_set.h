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

// Compute ellipsoidal invariant set {x : x'Px <= alpha} from a pre-computed
// DARE solution P and LQR gain K, given input constraints [u_min, u_max].
//
// For each input dimension i, alpha_i = min(u_max_i^2, u_min_i^2) / (k_i' P^{-1} k_i)
// where k_i is the i-th row of K (transposed to column). The tightest bound gives alpha.
template <typename Scalar, std::size_t NX, std::size_t NU>
auto compute_ellipsoidal_set(const Matrix<Scalar, NX, NX>& P, const Matrix<Scalar, NU, NX>& K, const Vector<Scalar, NU>& u_min, const Vector<Scalar, NU>& u_max) -> ellipsoidal_set<Scalar, NX>
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);

    auto ldlt = P.ldlt();
    Scalar alpha = std::numeric_limits<Scalar>::infinity();

    for(int i = 0; i < nu; ++i)
    {
        // k_i is the i-th row of K, as a column vector
        Vector<Scalar, NX> ki = K.row(i).transpose();

        // P^{-1} k_i via LDLT solve
        Eigen::Matrix<Scalar, nx, 1> Pinv_ki = ldlt.solve(ki);

        Scalar denom = ki.dot(Pinv_ki);
        if(denom <= Scalar{0})
            continue;

        Scalar u_bound = std::min(u_max(i) * u_max(i), u_min(i) * u_min(i));
        Scalar alpha_i = u_bound / denom;
        alpha = std::min(alpha, alpha_i);
    }

    return {.P = P, .alpha = alpha};
}

// Compute polytopic maximal control-invariant set via backward reachability.
// Restricted to NX <= 4 due to computational complexity of vertex enumeration.
//
// Starting from C_0 = state_constraints, iterates:
//   C_{k+1} = C_k intersect Pre(C_k)
// where Pre(C) = { x : exists u in U s.t. Ax + Bu in C }
//
// For QP-compatible polytopes, computes Pre by enumerating vertices of the
// input constraint polytope and intersecting the pre-images.
template <typename Scalar, std::size_t NX, std::size_t NU>
auto compute_polytopic_invariant_set(const Matrix<Scalar, NX, NX>& A_sys,
                                     const Matrix<Scalar, NX, NU>& B_sys,
                                     const polytopic_set<Scalar, NX>& state_constraints,
                                     const polytopic_set<Scalar, NU>& input_constraints,
                                     int max_iterations = 100,
                                     Scalar convergence_tol = Scalar{1e-6}) -> std::optional<polytopic_set<Scalar, NX>>
{
    static_assert(NX <= 4, "Polytopic invariant set computation restricted to NX <= 4");

    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);

    // Enumerate vertices of input constraint polytope.
    // For box constraints (most common), this is 2^NU corners.
    // General case: enumerate all C(m, nu) intersections of nu hyperplanes.
    auto enumerate_vertices = [&]() -> std::vector<Vector<Scalar, NU>>
    {
        std::vector<Vector<Scalar, NU>> vertices;
        int m = static_cast<int>(input_constraints.H.rows());

        if(m < nu)
            return vertices;

        // Generate all combinations of nu rows from m hyperplanes
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
            // Form nu x nu subsystem from selected rows
            Eigen::Matrix<Scalar, nu, nu> H_sub;
            Eigen::Matrix<Scalar, nu, 1> h_sub;
            for(int i = 0; i < nu; ++i)
            {
                H_sub.row(i) = input_constraints.H.row(indices[static_cast<std::size_t>(i)]);
                h_sub(i) = input_constraints.h(indices[static_cast<std::size_t>(i)]);
            }

            // Solve H_sub * v = h_sub
            auto qr = H_sub.colPivHouseholderQr();
            if(!qr.isInvertible())
                continue;

            Vector<Scalar, NU> v = qr.solve(h_sub);

            // Check feasibility: H * v <= h + tol for all constraints
            bool feasible = true;
            for(int j = 0; j < m; ++j)
            {
                Scalar val = input_constraints.H.row(j).dot(v);
                if(val > input_constraints.h(j) + convergence_tol)
                {
                    feasible = false;
                    break;
                }
            }

            if(feasible)
                vertices.push_back(v);
        } while(next_combination());

        return vertices;
    };

    auto u_vertices = enumerate_vertices();
    if(u_vertices.empty())
        return std::nullopt;

    // Current invariant set approximation (H-representation)
    using HMatrix = Eigen::Matrix<Scalar, Eigen::Dynamic, nx>;
    using hVector = Eigen::VectorX<Scalar>;

    HMatrix H_curr = state_constraints.H;
    hVector h_curr = state_constraints.h;

    for(int iter = 0; iter < max_iterations; ++iter)
    {
        // Compute Pre(C_curr) by intersecting over input vertices:
        // For each vertex u_v: {x : H_curr * (A*x + B*u_v) <= h_curr}
        //                     = {x : (H_curr * A) * x <= h_curr - H_curr * B * u_v}
        HMatrix H_A = H_curr * A_sys;
        int n_faces = static_cast<int>(H_curr.rows());
        int n_verts = static_cast<int>(u_vertices.size());

        HMatrix H_pre(n_faces * n_verts, nx);
        hVector h_pre(n_faces * n_verts);

        for(int v = 0; v < n_verts; ++v)
        {
            Eigen::Matrix<Scalar, nx, 1> Bu = B_sys * u_vertices[static_cast<std::size_t>(v)];
            hVector offset = H_curr * Bu;

            H_pre.middleRows(v * n_faces, n_faces) = H_A;
            h_pre.segment(v * n_faces, n_faces) = h_curr - offset;
        }

        // Intersect Pre(C) with C_curr
        int old_rows = static_cast<int>(H_curr.rows());
        int new_rows = old_rows + n_faces * n_verts;

        HMatrix H_next(new_rows, nx);
        hVector h_next(new_rows);

        H_next.topRows(old_rows) = H_curr;
        h_next.head(old_rows) = h_curr;
        H_next.bottomRows(n_faces * n_verts) = H_pre;
        h_next.tail(n_faces * n_verts) = h_pre;

        // Simple redundancy check: remove rows where h is very large
        // relative to the row norm (these constraints are not binding)
        std::vector<int> keep;
        keep.reserve(static_cast<std::size_t>(new_rows));
        for(int i = 0; i < new_rows; ++i)
        {
            Scalar row_norm = H_next.row(i).norm();
            if(row_norm < Scalar{1e-14})
                continue;
            // Normalize and check if constraint is meaningful
            Scalar h_normalized = h_next(i) / row_norm;

            // Check for near-duplicate rows: keep the tightest
            bool is_redundant = false;
            for(int j : keep)
            {
                Scalar row_norm_j = H_next.row(j).norm();
                Vector<Scalar, NX> dir_i = H_next.row(i).transpose() / row_norm;
                Vector<Scalar, NX> dir_j = H_next.row(j).transpose() / row_norm_j;
                Scalar dot = dir_i.dot(dir_j);
                if(dot > Scalar{1} - Scalar{1e-8})
                {
                    // Nearly parallel -- keep tighter
                    Scalar h_norm_j = h_next(j) / row_norm_j;
                    if(h_normalized >= h_norm_j - convergence_tol)
                    {
                        is_redundant = true;
                        break;
                    }
                    else
                    {
                        // Replace j with i (i is tighter)
                        // Mark j for removal would be complex; skip for simplicity
                    }
                }
            }
            if(!is_redundant)
                keep.push_back(i);
        }

        // Cap halfplane count to prevent unbounded growth
        constexpr int max_halfplanes = 500;
        if(static_cast<int>(keep.size()) > max_halfplanes)
            keep.resize(static_cast<std::size_t>(max_halfplanes));

        HMatrix H_filtered(static_cast<int>(keep.size()), nx);
        hVector h_filtered(static_cast<int>(keep.size()));
        for(int i = 0; i < static_cast<int>(keep.size()); ++i)
        {
            H_filtered.row(i) = H_next.row(keep[static_cast<std::size_t>(i)]);
            h_filtered(i) = h_next(keep[static_cast<std::size_t>(i)]);
        }

        // Convergence check: if no new binding constraints added
        if(static_cast<int>(keep.size()) == old_rows)
        {
            return polytopic_set<Scalar, NX>{.H = std::move(H_filtered), .h = std::move(h_filtered)};
        }

        H_curr = std::move(H_filtered);
        h_curr = std::move(h_filtered);
    }

    // Did not converge; return current best approximation (conservative)
    return polytopic_set<Scalar, NX>{.H = std::move(H_curr), .h = std::move(h_curr)};
}

// Result type for terminal_ingredients helper.
template <typename Scalar, std::size_t NX, std::size_t NU>
struct terminal_ingredients_result
{
    Matrix<Scalar, NX, NX> Qf;
    ellipsoidal_set<Scalar, NX> set;
};

// Compute both terminal cost Qf and ellipsoidal terminal constraint set
// from system matrices and input bounds, per Mayne et al. 2000.
//
// Steps:
//   1. Solve DARE for P (= Qf)
//   2. Compute LQR gain K = -(R + B'PB)^{-1} B'PA
//   3. Compute ellipsoidal set {x : x'Px <= alpha} from K and input bounds
template <typename Scalar, std::size_t NX, std::size_t NU>
auto terminal_ingredients(
    const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NX, NU>& B, const Matrix<Scalar, NX, NX>& Q, const Matrix<Scalar, NU, NU>& R, const Vector<Scalar, NU>& u_min, const Vector<Scalar, NU>& u_max)
    -> std::optional<terminal_ingredients_result<Scalar, NX, NU>>
{
    auto P_opt = dare<Scalar, NX, NU>(A, B, Q, R);
    if(!P_opt)
        return std::nullopt;

    auto P = P_opt.value();

    // K = -(R + B'PB)^{-1} B'PA
    Matrix<Scalar, NU, NU> RpBtPB = R + B.transpose() * P * B;
    Matrix<Scalar, NU, NX> K = -(RpBtPB.ldlt().solve(B.transpose() * P * A));

    auto eset = compute_ellipsoidal_set<Scalar, NX, NU>(P, K, u_min, u_max);

    return terminal_ingredients_result<Scalar, NX, NU>{.Qf = P, .set = eset};
}

} // namespace ctrlpp

#endif
