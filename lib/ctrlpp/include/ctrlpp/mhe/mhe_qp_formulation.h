#ifndef HPP_GUARD_CTRLPP_MHE_MHE_QP_FORMULATION_H
#define HPP_GUARD_CTRLPP_MHE_MHE_QP_FORMULATION_H

#include "ctrlpp/types.h"

#include "ctrlpp/mpc/qp_types.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <array>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace ctrlpp::detail {

/// Compute dimensions for the MHE QP.
///
/// Decision vector layout: z = [x_0, x_1, ..., x_N, s_0, ..., s_N]
/// where s_k are slack variables for soft box constraints (present only when
/// soft_constraints=true and bounds exist).
///
/// Constraint layout:
///   Block 1: Dynamics equality  -- N rows of NX (x_{k+1} = A_k x_k + B_k u_k)
///   Block 2: State box bounds   -- (N+1)*NX rows (identity on each x_k, with -I on slack)
///   Block 3: Residual bounds    -- (N+1)*NY rows per direction (when residual_bound present)
struct mhe_qp_dims {
    int n_states;
    int n_slack;
    int n_dec;
    int n_dyn;
    int n_box;
    int n_residual;
    int n_con;
};

template<std::size_t NX, std::size_t NY>
[[nodiscard]] constexpr auto compute_mhe_dims(
    std::size_t N,
    bool has_box_bounds,
    bool has_soft_constraints,
    bool has_residual_bounds) -> mhe_qp_dims
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int ny = static_cast<int>(NY);
    int Ni = static_cast<int>(N);

    int n_states = (Ni + 1) * nx;
    bool need_slack = has_soft_constraints && has_box_bounds;
    int n_slack = need_slack ? (Ni + 1) * nx : 0;
    int n_dec = n_states + n_slack;

    int n_dyn = Ni * nx;
    int n_box = has_box_bounds ? (Ni + 1) * nx : 0;
    int n_residual = has_residual_bounds ? (Ni + 1) * ny : 0;
    int n_con = n_dyn + n_box + n_residual;

    return {n_states, n_slack, n_dec, n_dyn, n_box, n_residual, n_con};
}

/// Build the constant parts of the MHE QP: Hessian P and constraint matrix A.
///
/// The Hessian encodes:
///   - Arrival cost: arrival_weight * P_arr_inv on x_0
///   - Process noise: Q_inv on each (x_{k+1} - A_k x_k - B_k u_k) expanded as
///     quadratic form over x_k and x_{k+1} with cross-terms
///   - Measurement: H_k^T R_inv H_k added to diagonal blocks
///   - Slack regularization: small diagonal for slack variables
///
/// The constraint matrix encodes:
///   - Dynamics equality: -A_k x_k + I x_{k+1} = B_k u_k (RHS in bounds)
///   - Box constraints: I x_k - I s_k (bounds = x_min, x_max)
///   - Residual bounds: H_k x_k (bounds = z_k +/- threshold)
template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
[[nodiscard]] auto build_mhe_qp_structure(
    std::size_t N,
    Scalar arrival_weight,
    const Matrix<Scalar, NX, NX>& P_arr_inv,
    const Matrix<Scalar, NX, NX>& Q_inv,
    const Matrix<Scalar, NY, NY>& R_inv,
    const std::array<Matrix<Scalar, NX, NX>, 1>& A_lin,
    const std::array<Matrix<Scalar, NY, NX>, 1>& H_lin,
    bool has_box_bounds,
    bool has_soft_constraints,
    Scalar soft_penalty,
    bool has_residual_bounds) -> qp_problem<Scalar>
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int ny = static_cast<int>(NY);
    int Ni = static_cast<int>(N);

    auto dims = compute_mhe_dims<NX, NY>(N, has_box_bounds, has_soft_constraints, has_residual_bounds);

    // --- Build Hessian P (upper triangular) ---
    std::vector<Eigen::Triplet<Scalar>> p_trips;
    p_trips.reserve(static_cast<std::size_t>(
        (Ni + 1) * nx * (nx + 1) / 2 * 3 + dims.n_slack));

    const auto& A = A_lin[0];
    const auto& H = H_lin[0];

    // Pre-compute recurring blocks
    Matrix<Scalar, NX, NX> AtQinvA = (A.transpose() * Q_inv * A).eval();
    Matrix<Scalar, NX, NX> AtQinv  = (A.transpose() * Q_inv).eval();
    Matrix<Scalar, NX, NX> HtRinvH = (H.transpose() * R_inv * H).eval();

    auto add_upper_block = [&](int row_off, int col_off, const auto& block, int rows, int cols) {
        for (int i = 0; i < rows; ++i) {
            int start_j = (row_off == col_off) ? i : 0;
            for (int j = start_j; j < cols; ++j) {
                if (block(i, j) != Scalar{0})
                    p_trips.emplace_back(row_off + i, col_off + j, block(i, j));
            }
        }
    };

    // Block (0,0): arrival cost + Q_inv from process noise + measurement
    // For k=0: arrival_weight * P_arr_inv + A^T Q_inv A + H^T R_inv H
    {
        Matrix<Scalar, NX, NX> blk = arrival_weight * P_arr_inv + AtQinvA + HtRinvH;
        add_upper_block(0, 0, blk, nx, nx);
    }

    // Blocks (k,k) for k=1..N-1: Q_inv + A^T Q_inv A + H^T R_inv H
    // (receives process noise from both step k-1->k and k->k+1)
    for (int k = 1; k < Ni; ++k) {
        Matrix<Scalar, NX, NX> blk = Q_inv + AtQinvA + HtRinvH;
        int off = k * nx;
        add_upper_block(off, off, blk, nx, nx);
    }

    // Block (N,N): Q_inv + H^T R_inv H (only receives process noise from step N-1->N)
    {
        Matrix<Scalar, NX, NX> blk = Q_inv + HtRinvH;
        int off = Ni * nx;
        add_upper_block(off, off, blk, nx, nx);
    }

    // Cross-terms (k, k+1) for k=0..N-1: -A^T Q_inv (upper triangular part of full Hessian)
    Matrix<Scalar, NX, NX> cross = (-AtQinv).eval();
    for (int k = 0; k < Ni; ++k) {
        int row_off = k * nx;
        int col_off = (k + 1) * nx;
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < nx; ++j) {
                if (cross(i, j) != Scalar{0})
                    p_trips.emplace_back(row_off + i, col_off + j, cross(i, j));
            }
        }
    }

    // Slack regularization (small positive diagonal)
    if (dims.n_slack > 0) {
        constexpr Scalar slack_reg{1e-6};
        int slack_off = dims.n_states;
        for (int i = 0; i < dims.n_slack; ++i) {
            p_trips.emplace_back(slack_off + i, slack_off + i, slack_reg);
        }
    }

    Eigen::SparseMatrix<Scalar, Eigen::ColMajor> P(dims.n_dec, dims.n_dec);
    P.setFromTriplets(p_trips.begin(), p_trips.end());

    // --- Build constraint matrix A_con ---
    std::vector<Eigen::Triplet<Scalar>> a_trips;
    a_trips.reserve(static_cast<std::size_t>(
        Ni * (nx * nx + nx) + dims.n_box * 2 + dims.n_residual * nx));

    int row = 0;

    // Block 1: Dynamics equality -- N rows of NX each
    // -A x_k + I x_{k+1} = B u_k (RHS set in bounds)
    for (int k = 0; k < Ni; ++k) {
        int xk_off = k * nx;
        int xk1_off = (k + 1) * nx;

        // -A on x_k
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < nx; ++j) {
                if (A(i, j) != Scalar{0})
                    a_trips.emplace_back(row + i, xk_off + j, -A(i, j));
            }
        }
        // I on x_{k+1}
        for (int i = 0; i < nx; ++i)
            a_trips.emplace_back(row + i, xk1_off + i, Scalar{1});

        row += nx;
    }

    // Block 2: Box constraints on x_0 .. x_N
    if (has_box_bounds) {
        bool need_slack = has_soft_constraints;
        for (int k = 0; k <= Ni; ++k) {
            int xk_off = k * nx;
            for (int i = 0; i < nx; ++i)
                a_trips.emplace_back(row + i, xk_off + i, Scalar{1});

            if (need_slack) {
                int sk_off = dims.n_states + k * nx;
                for (int i = 0; i < nx; ++i)
                    a_trips.emplace_back(row + i, sk_off + i, Scalar{-1});
            }
            row += nx;
        }
    }

    // Block 3: Residual bounds -- H_k x_k (bounds set per-solve)
    if (has_residual_bounds) {
        for (int k = 0; k <= Ni; ++k) {
            int xk_off = k * nx;
            for (int i = 0; i < ny; ++i) {
                for (int j = 0; j < nx; ++j) {
                    if (H(i, j) != Scalar{0})
                        a_trips.emplace_back(row + i, xk_off + j, H(i, j));
                }
            }
            row += ny;
        }
    }

    Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A_con(dims.n_con, dims.n_dec);
    A_con.setFromTriplets(a_trips.begin(), a_trips.end());

    // --- Build initial q, l, u vectors (placeholder, updated each solve) ---
    Eigen::VectorX<Scalar> q = Eigen::VectorX<Scalar>::Zero(dims.n_dec);

    // Slack linear cost (L1 penalty)
    if (dims.n_slack > 0) {
        int slack_off = dims.n_states;
        for (int i = 0; i < dims.n_slack; ++i)
            q(slack_off + i) = soft_penalty;
    }

    Eigen::VectorX<Scalar> l = Eigen::VectorX<Scalar>::Zero(dims.n_con);
    Eigen::VectorX<Scalar> u = Eigen::VectorX<Scalar>::Zero(dims.n_con);

    return {std::move(P), std::move(q), std::move(A_con), std::move(l), std::move(u)};
}

/// Build per-solve QP update vectors (q, l, u) from current window data.
///
/// Called each time the MHE window shifts and a new solve is needed.
template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
[[nodiscard]] auto build_mhe_qp_update(
    std::size_t N,
    Scalar arrival_weight,
    const Matrix<Scalar, NX, NX>& P_arr_inv,
    const Matrix<Scalar, NX, NX>& Q_inv,
    const Matrix<Scalar, NY, NY>& R_inv,
    const Matrix<Scalar, NX, NX>& A_lin,
    const Matrix<Scalar, NX, NU>& B_lin,
    const Matrix<Scalar, NY, NX>& H_lin,
    const Vector<Scalar, NX>& x_arrival,
    std::span<const Vector<Scalar, NU>> u_buf,
    std::span<const Vector<Scalar, NY>> z_buf,
    bool has_box_bounds,
    bool has_soft_constraints,
    Scalar soft_penalty,
    const std::optional<Vector<Scalar, NX>>& x_min,
    const std::optional<Vector<Scalar, NX>>& x_max,
    bool has_residual_bounds,
    const std::optional<Vector<Scalar, NY>>& residual_bound,
    const Eigen::VectorX<Scalar>& warm_x,
    const Eigen::VectorX<Scalar>& warm_y) -> qp_update<Scalar>
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);
    constexpr int ny = static_cast<int>(NY);
    constexpr auto inf = std::numeric_limits<Scalar>::infinity();
    int Ni = static_cast<int>(N);

    auto dims = compute_mhe_dims<NX, NY>(N, has_box_bounds, has_soft_constraints, has_residual_bounds);

    // Pre-compute
    Matrix<Scalar, NX, NX> AtQinv = (A_lin.transpose() * Q_inv).eval();
    Matrix<Scalar, NX, NY> HtRinv = (H_lin.transpose() * R_inv).eval();

    // --- Build linear cost q ---
    Eigen::VectorX<Scalar> q = Eigen::VectorX<Scalar>::Zero(dims.n_dec);

    // Arrival cost linear term on x_0: -arrival_weight * P_arr_inv * x_arrival
    q.template segment<nx>(0) = -arrival_weight * P_arr_inv * x_arrival;

    // Process noise cross-terms from known inputs: A^T Q_inv B u_k
    // For x_k (k=0..N-1): add A^T Q_inv B u_k
    // For x_{k+1} (k=0..N-1): add -Q_inv B u_k
    for (int k = 0; k < Ni; ++k) {
        Vector<Scalar, NX> Buk = B_lin * u_buf[static_cast<std::size_t>(k)];
        Vector<Scalar, NX> QinvBuk = (Q_inv * Buk).eval();

        q.segment(k * nx, nx) += AtQinv * Buk;
        q.segment((k + 1) * nx, nx) -= QinvBuk;
    }

    // Measurement terms: -H^T R_inv z_k on each x_k
    for (int k = 0; k <= Ni; ++k) {
        q.segment(k * nx, nx) -= HtRinv * z_buf[static_cast<std::size_t>(k)];
    }

    // Slack L1 penalty (linear cost)
    if (dims.n_slack > 0) {
        int slack_off = dims.n_states;
        for (int i = 0; i < dims.n_slack; ++i)
            q(slack_off + i) = soft_penalty;
    }

    // --- Build bounds l, u ---
    Eigen::VectorX<Scalar> l(dims.n_con);
    Eigen::VectorX<Scalar> u(dims.n_con);

    int row = 0;

    // Dynamics equality: -A x_k + x_{k+1} = B u_k
    for (int k = 0; k < Ni; ++k) {
        Vector<Scalar, NX> rhs = B_lin * u_buf[static_cast<std::size_t>(k)];
        l.segment(row, nx) = rhs;
        u.segment(row, nx) = rhs;
        row += nx;
    }

    // Box constraints on x_0 .. x_N
    if (has_box_bounds) {
        Vector<Scalar, NX> lb = x_min.has_value() ? x_min.value()
            : Vector<Scalar, NX>::Constant(-inf);
        Vector<Scalar, NX> ub = x_max.has_value() ? x_max.value()
            : Vector<Scalar, NX>::Constant(inf);
        for (int k = 0; k <= Ni; ++k) {
            l.segment(row, nx) = lb;
            u.segment(row, nx) = ub;
            row += nx;
        }
    }

    // Residual bounds: z_k,i - threshold_i <= H_i x_k <= z_k,i + threshold_i
    if (has_residual_bounds && residual_bound.has_value()) {
        const auto& thresh = residual_bound.value();
        for (int k = 0; k <= Ni; ++k) {
            const auto& zk = z_buf[static_cast<std::size_t>(k)];
            for (int i = 0; i < ny; ++i) {
                l(row + i) = zk(i) - thresh(i);
                u(row + i) = zk(i) + thresh(i);
            }
            row += ny;
        }
    }

    return {std::move(q), std::move(l), std::move(u), warm_x, warm_y};
}

}

#endif
