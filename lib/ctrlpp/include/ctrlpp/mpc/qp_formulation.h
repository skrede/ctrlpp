#ifndef HPP_GUARD_CTRLPP_MPC_QP_FORMULATION_H
#define HPP_GUARD_CTRLPP_MPC_QP_FORMULATION_H

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/types.h"

#include <Eigen/Sparse>

#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace ctrlpp::detail {

// Build cost Hessian P (upper triangular, sparse).
// P = blkdiag(Q, ..., Q, Qf, R, ..., R) for states and inputs.
// With soft constraints: append blkdiag(rho*I, ..., rho*I) for slack variables.
template<typename Scalar, std::size_t NX, std::size_t NU>
[[nodiscard]] auto build_cost_matrix(
    int N,
    const Matrix<Scalar, NX, NX>& Q,
    const Matrix<Scalar, NU, NU>& R,
    const Matrix<Scalar, NX, NX>& Qf,
    bool has_soft_constraints,
    Scalar soft_penalty) -> Eigen::SparseMatrix<Scalar, Eigen::ColMajor>
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);

    int n_x_total = (N + 1) * nx;
    int n_u_total = N * nu;
    int n_slack = has_soft_constraints ? N * nx : 0;
    int n_dec = n_x_total + n_u_total + n_slack;

    std::vector<Eigen::Triplet<Scalar>> triplets;
    triplets.reserve(static_cast<std::size_t>(
        (N + 1) * nx * (nx + 1) / 2 + N * nu * (nu + 1) / 2 + n_slack));

    // Q blocks for x_0 .. x_{N-1}
    for (int k = 0; k < N; ++k) {
        int offset = k * nx;
        for (int i = 0; i < nx; ++i) {
            for (int j = i; j < nx; ++j) {
                if (Q(i, j) != Scalar{0})
                    triplets.emplace_back(offset + i, offset + j, Q(i, j));
            }
        }
    }

    // Qf block for x_N
    {
        int offset = N * nx;
        for (int i = 0; i < nx; ++i) {
            for (int j = i; j < nx; ++j) {
                if (Qf(i, j) != Scalar{0})
                    triplets.emplace_back(offset + i, offset + j, Qf(i, j));
            }
        }
    }

    // R blocks for u_0 .. u_{N-1}
    for (int k = 0; k < N; ++k) {
        int offset = n_x_total + k * nu;
        for (int i = 0; i < nu; ++i) {
            for (int j = i; j < nu; ++j) {
                if (R(i, j) != Scalar{0})
                    triplets.emplace_back(offset + i, offset + j, R(i, j));
            }
        }
    }

    // Slack penalty blocks: rho * I for each of N steps
    if (has_soft_constraints) {
        int slack_offset = n_x_total + n_u_total;
        for (int k = 0; k < N; ++k) {
            for (int i = 0; i < nx; ++i) {
                triplets.emplace_back(
                    slack_offset + k * nx + i,
                    slack_offset + k * nx + i,
                    soft_penalty);
            }
        }
    }

    Eigen::SparseMatrix<Scalar, Eigen::ColMajor> P(n_dec, n_dec);
    P.setFromTriplets(triplets.begin(), triplets.end());
    return P;
}

// Build constraint matrix A (sparse).
// Row blocks:
//   1. Dynamics equality: (N+1)*NX rows
//   2. State bounds: N*NX rows (x_1 .. x_N, skip x_0)
//   3. Input bounds: N*NU rows
//   4. Rate bounds: N*NU rows (if has_rate_bounds)
template<typename Scalar, std::size_t NX, std::size_t NU>
[[nodiscard]] auto build_constraint_matrix(
    int N,
    const Matrix<Scalar, NX, NX>& A_sys,
    const Matrix<Scalar, NX, NU>& B_sys,
    bool has_state_bounds,
    bool has_soft_constraints,
    bool has_input_bounds,
    bool has_rate_bounds) -> Eigen::SparseMatrix<Scalar, Eigen::ColMajor>
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);

    int n_x_total = (N + 1) * nx;
    int n_u_total = N * nu;
    int n_slack = has_soft_constraints ? N * nx : 0;
    int n_dec = n_x_total + n_u_total + n_slack;

    int n_dyn = (N + 1) * nx;
    int n_state = has_state_bounds ? N * nx : 0;
    int n_input = has_input_bounds ? N * nu : 0;
    int n_rate = has_rate_bounds ? N * nu : 0;
    int n_con = n_dyn + n_state + n_input + n_rate;

    std::vector<Eigen::Triplet<Scalar>> triplets;
    triplets.reserve(static_cast<std::size_t>(
        (N + 1) * nx + N * nx * nx + N * nx * nu
        + n_state + n_slack + n_input + 2 * n_rate));

    int row = 0;

    // Block 1: Dynamics equality constraints
    // Row 0: I * x_0 = x0  (identity on x_0)
    for (int i = 0; i < nx; ++i)
        triplets.emplace_back(row + i, i, Scalar{1});
    row += nx;

    // Rows k=1..N: -A_sys * x_{k-1} + I * x_k - B_sys * u_{k-1} = 0
    for (int k = 1; k <= N; ++k) {
        int x_k_offset = k * nx;
        int x_km1_offset = (k - 1) * nx;
        int u_km1_offset = n_x_total + (k - 1) * nu;

        // -A_sys on x_{k-1}
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < nx; ++j) {
                if (A_sys(i, j) != Scalar{0})
                    triplets.emplace_back(row + i, x_km1_offset + j, -A_sys(i, j));
            }
        }

        // I on x_k
        for (int i = 0; i < nx; ++i)
            triplets.emplace_back(row + i, x_k_offset + i, Scalar{1});

        // -B_sys on u_{k-1}
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < nu; ++j) {
                if (B_sys(i, j) != Scalar{0})
                    triplets.emplace_back(row + i, u_km1_offset + j, -B_sys(i, j));
            }
        }

        row += nx;
    }

    // Block 2: State bounds on x_1 .. x_N
    if (has_state_bounds) {
        for (int k = 1; k <= N; ++k) {
            int x_k_offset = k * nx;
            for (int i = 0; i < nx; ++i) {
                triplets.emplace_back(row + i, x_k_offset + i, Scalar{1});
            }
            // Soft constraints: -I on corresponding slack column
            if (has_soft_constraints) {
                int slack_offset = n_x_total + n_u_total + (k - 1) * nx;
                for (int i = 0; i < nx; ++i)
                    triplets.emplace_back(row + i, slack_offset + i, Scalar{-1});
            }
            row += nx;
        }
    }

    // Block 3: Input bounds on u_0 .. u_{N-1}
    if (has_input_bounds) {
        for (int k = 0; k < N; ++k) {
            int u_k_offset = n_x_total + k * nu;
            for (int i = 0; i < nu; ++i)
                triplets.emplace_back(row + i, u_k_offset + i, Scalar{1});
            row += nu;
        }
    }

    // Block 4: Rate bounds (delta-u constraints)
    if (has_rate_bounds) {
        for (int k = 0; k < N; ++k) {
            int u_k_offset = n_x_total + k * nu;
            // u_k
            for (int i = 0; i < nu; ++i)
                triplets.emplace_back(row + i, u_k_offset + i, Scalar{1});
            // -u_{k-1} (for k > 0)
            if (k > 0) {
                int u_km1_offset = n_x_total + (k - 1) * nu;
                for (int i = 0; i < nu; ++i)
                    triplets.emplace_back(row + i, u_km1_offset + i, Scalar{-1});
            }
            row += nu;
        }
    }

    Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A(n_con, n_dec);
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A;
}

// Build initial bound vectors l, u for the QP constraints.
template<typename Scalar, std::size_t NX, std::size_t NU>
[[nodiscard]] auto build_bounds_vectors(
    int N,
    const Vector<Scalar, NX>& x0,
    const std::optional<Vector<Scalar, NX>>& x_min,
    const std::optional<Vector<Scalar, NX>>& x_max,
    const std::optional<Vector<Scalar, NU>>& u_min,
    const std::optional<Vector<Scalar, NU>>& u_max,
    const std::optional<Vector<Scalar, NU>>& du_max,
    bool has_soft_constraints) -> std::pair<Eigen::VectorX<Scalar>, Eigen::VectorX<Scalar>>
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);
    constexpr auto inf = std::numeric_limits<Scalar>::infinity();

    bool has_state_bounds = x_min.has_value() || x_max.has_value();
    bool has_input_bounds = u_min.has_value() || u_max.has_value();
    bool has_rate_bounds = du_max.has_value();

    int n_dyn = (N + 1) * nx;
    int n_state = has_state_bounds ? N * nx : 0;
    int n_input = has_input_bounds ? N * nu : 0;
    int n_rate = has_rate_bounds ? N * nu : 0;
    int n_con = n_dyn + n_state + n_input + n_rate;

    Eigen::VectorX<Scalar> l(n_con);
    Eigen::VectorX<Scalar> u(n_con);

    int row = 0;

    // Dynamics equality: l = u = [x0, 0, 0, ...]
    l.segment(row, nx) = x0;
    u.segment(row, nx) = x0;
    row += nx;
    for (int k = 1; k <= N; ++k) {
        l.segment(row, nx).setZero();
        u.segment(row, nx).setZero();
        row += nx;
    }

    // State bounds on x_1 .. x_N
    // With soft constraints: x_min - s <= x <= x_max + s becomes
    // bounds on (x - s): x_min <= x - s, x - s <= x_max
    // (the slack is subtracted in the constraint matrix, so bounds stay as x_min/x_max)
    if (has_state_bounds) {
        Vector<Scalar, NX> lb = x_min.has_value() ? x_min.value()
            : Vector<Scalar, NX>::Constant(-inf);
        Vector<Scalar, NX> ub = x_max.has_value() ? x_max.value()
            : Vector<Scalar, NX>::Constant(inf);
        for (int k = 0; k < N; ++k) {
            l.segment(row, nx) = lb;
            u.segment(row, nx) = ub;
            row += nx;
        }
    }

    // Input bounds
    if (has_input_bounds) {
        Vector<Scalar, NU> lb = u_min.has_value() ? u_min.value()
            : Vector<Scalar, NU>::Constant(-inf);
        Vector<Scalar, NU> ub = u_max.has_value() ? u_max.value()
            : Vector<Scalar, NU>::Constant(inf);
        for (int k = 0; k < N; ++k) {
            l.segment(row, nu) = lb;
            u.segment(row, nu) = ub;
            row += nu;
        }
    }

    // Rate bounds
    if (has_rate_bounds) {
        for (int k = 0; k < N; ++k) {
            l.segment(row, nu) = -du_max.value();
            u.segment(row, nu) = du_max.value();
            row += nu;
        }
    }

    return {l, u};
}

// Build cost vector q from reference trajectory.
// q = [-Q*xref_0, ..., -Q*xref_{N-1}, -Qf*xref_N, 0, ..., 0]
// For origin regulation (no reference), q is all zeros.
template<typename Scalar, std::size_t NX, std::size_t NU>
[[nodiscard]] auto build_cost_vector(
    int N,
    int n_dec,
    const Matrix<Scalar, NX, NX>& Q,
    const Matrix<Scalar, NX, NX>& Qf) -> Eigen::VectorX<Scalar>
{
    // Origin regulation: all zeros
    return Eigen::VectorX<Scalar>::Zero(n_dec);
}

// Overload: single setpoint reference broadcast to all steps.
template<typename Scalar, std::size_t NX, std::size_t NU>
[[nodiscard]] auto build_cost_vector(
    int N,
    int n_dec,
    const Matrix<Scalar, NX, NX>& Q,
    const Matrix<Scalar, NX, NX>& Qf,
    const Vector<Scalar, NX>& x_ref) -> Eigen::VectorX<Scalar>
{
    constexpr int nx = static_cast<int>(NX);

    Eigen::VectorX<Scalar> q = Eigen::VectorX<Scalar>::Zero(n_dec);

    // -Q * x_ref for steps 0..N-1
    Vector<Scalar, NX> Qxr = Q * x_ref;
    for (int k = 0; k < N; ++k)
        q.segment(k * nx, nx) = -Qxr;

    // -Qf * x_ref for step N
    q.segment(N * nx, nx) = -(Qf * x_ref);

    return q;
}

// Overload: per-step reference trajectory (span of size N+1).
template<typename Scalar, std::size_t NX, std::size_t NU>
[[nodiscard]] auto build_cost_vector(
    int N,
    int n_dec,
    const Matrix<Scalar, NX, NX>& Q,
    const Matrix<Scalar, NX, NX>& Qf,
    std::span<const Vector<Scalar, NX>> x_ref) -> Eigen::VectorX<Scalar>
{
    constexpr int nx = static_cast<int>(NX);

    Eigen::VectorX<Scalar> q = Eigen::VectorX<Scalar>::Zero(n_dec);

    // -Q * x_ref_k for steps 0..N-1
    for (int k = 0; k < N; ++k)
        q.segment(k * nx, nx) = -(Q * x_ref[static_cast<std::size_t>(k)]);

    // -Qf * x_ref_N for step N
    q.segment(N * nx, nx) = -(Qf * x_ref[static_cast<std::size_t>(N)]);

    return q;
}

}

#endif
