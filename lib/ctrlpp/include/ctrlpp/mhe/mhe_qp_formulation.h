#ifndef HPP_GUARD_CTRLPP_MHE_MHE_QP_FORMULATION_H
#define HPP_GUARD_CTRLPP_MHE_MHE_QP_FORMULATION_H

/// @brief QP formulation for linear Moving Horizon Estimation.
///
/// @cite rao2003 -- Rao et al., "Constrained State Estimation for Nonlinear Discrete-Time Systems", 2003

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

namespace ctrlpp::detail
{

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
struct mhe_qp_dims
{
    int n_states;
    int n_slack;
    int n_dec;
    int n_dyn;
    int n_box;
    int n_residual;
    int n_con;
};

template <std::size_t NX, std::size_t NY>
[[nodiscard]] constexpr auto compute_mhe_dims(std::size_t N, bool has_box_bounds, bool has_soft_constraints, bool has_residual_bounds) -> mhe_qp_dims
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

/// Add upper-triangular entries from a dense block to the triplet list.
template <typename Scalar>
inline void add_upper_block(std::vector<Eigen::Triplet<Scalar>>& trips,
                            int row_off, int col_off, const auto& block, int rows, int cols)
{
    for(int i = 0; i < rows; ++i)
    {
        int start_j = (row_off == col_off) ? i : 0;
        for(int j = start_j; j < cols; ++j)
            if(block(i, j) != Scalar{0})
                trips.emplace_back(row_off + i, col_off + j, block(i, j));
    }
}

/// @cite rao2003 -- Hessian diagonal blocks from arrival cost, process noise, and measurement noise
template <typename Scalar, std::size_t NX>
inline void build_mhe_hessian_diagonal(std::vector<Eigen::Triplet<Scalar>>& trips,
                                       int Ni,
                                       Scalar arrival_weight,
                                       const Matrix<Scalar, NX, NX>& P_arr_inv,
                                       const Matrix<Scalar, NX, NX>& AtQinvA,
                                       const Matrix<Scalar, NX, NX>& HtRinvH,
                                       const Matrix<Scalar, NX, NX>& Q_inv)
{
    constexpr int nx = static_cast<int>(NX);

    // Block (0,0): arrival + process noise + measurement
    Matrix<Scalar, NX, NX> blk0 = arrival_weight * P_arr_inv + AtQinvA + HtRinvH;
    add_upper_block<Scalar>(trips, 0, 0, blk0, nx, nx);

    // Blocks (k,k) for k=1..N-1: receives process noise from both sides
    Matrix<Scalar, NX, NX> blk_mid = Q_inv + AtQinvA + HtRinvH;
    for(int k = 1; k < Ni; ++k)
        add_upper_block<Scalar>(trips, k * nx, k * nx, blk_mid, nx, nx);

    // Block (N,N): only receives process noise from step N-1->N
    Matrix<Scalar, NX, NX> blk_last = Q_inv + HtRinvH;
    add_upper_block<Scalar>(trips, Ni * nx, Ni * nx, blk_last, nx, nx);
}

/// @cite rao2003 -- Hessian cross-term blocks -A^T Q^{-1} between consecutive states
template <typename Scalar, std::size_t NX>
inline void build_mhe_hessian_cross_terms(std::vector<Eigen::Triplet<Scalar>>& trips,
                                          int Ni,
                                          const Matrix<Scalar, NX, NX>& AtQinv)
{
    constexpr int nx = static_cast<int>(NX);
    Matrix<Scalar, NX, NX> cross = (-AtQinv).eval();
    for(int k = 0; k < Ni; ++k)
    {
        int row_off = k * nx;
        int col_off = (k + 1) * nx;
        for(int i = 0; i < nx; ++i)
            for(int j = 0; j < nx; ++j)
                if(cross(i, j) != Scalar{0})
                    trips.emplace_back(row_off + i, col_off + j, cross(i, j));
    }
}

/// Build MHE Hessian matrix from cost components.
template <typename Scalar, std::size_t NX, std::size_t NY>
[[nodiscard]] auto build_mhe_hessian(const mhe_qp_dims& dims,
                                     int Ni,
                                     Scalar arrival_weight,
                                     const Matrix<Scalar, NX, NX>& P_arr_inv,
                                     const Matrix<Scalar, NX, NX>& Q_inv,
                                     const Matrix<Scalar, NY, NY>& R_inv,
                                     const Matrix<Scalar, NX, NX>& A,
                                     const Matrix<Scalar, NY, NX>& H) -> Eigen::SparseMatrix<Scalar, Eigen::ColMajor>
{
    constexpr int nx = static_cast<int>(NX);

    std::vector<Eigen::Triplet<Scalar>> p_trips;
    p_trips.reserve(static_cast<std::size_t>((Ni + 1) * nx * (nx + 1) / 2 * 3 + dims.n_slack));

    Matrix<Scalar, NX, NX> AtQinvA = (A.transpose() * Q_inv * A).eval();
    Matrix<Scalar, NX, NX> AtQinv = (A.transpose() * Q_inv).eval();
    Matrix<Scalar, NX, NX> HtRinvH = (H.transpose() * R_inv * H).eval();

    build_mhe_hessian_diagonal<Scalar, NX>(p_trips, Ni, arrival_weight, P_arr_inv, AtQinvA, HtRinvH, Q_inv);
    build_mhe_hessian_cross_terms<Scalar, NX>(p_trips, Ni, AtQinv);

    if(dims.n_slack > 0)
    {
        constexpr Scalar slack_reg{1e-6};
        int slack_off = dims.n_states;
        for(int i = 0; i < dims.n_slack; ++i)
            p_trips.emplace_back(slack_off + i, slack_off + i, slack_reg);
    }

    Eigen::SparseMatrix<Scalar, Eigen::ColMajor> P(dims.n_dec, dims.n_dec);
    P.setFromTriplets(p_trips.begin(), p_trips.end());
    return P;
}

/// Build MHE dynamics constraint triplets: -A x_k + I x_{k+1} = B u_k
template <typename Scalar, std::size_t NX>
inline void add_mhe_dynamics_triplets(std::vector<Eigen::Triplet<Scalar>>& trips, int& row, int Ni, const Matrix<Scalar, NX, NX>& A)
{
    constexpr int nx = static_cast<int>(NX);
    for(int k = 0; k < Ni; ++k)
    {
        int xk_off = k * nx;
        int xk1_off = (k + 1) * nx;

        for(int i = 0; i < nx; ++i)
            for(int j = 0; j < nx; ++j)
                if(A(i, j) != Scalar{0})
                    trips.emplace_back(row + i, xk_off + j, -A(i, j));
        for(int i = 0; i < nx; ++i)
            trips.emplace_back(row + i, xk1_off + i, Scalar{1});
        row += nx;
    }
}

/// Build MHE box constraint triplets: I x_k - I s_k
template <typename Scalar, std::size_t NX>
inline void add_mhe_box_triplets(std::vector<Eigen::Triplet<Scalar>>& trips, int& row, int Ni, int n_states, bool need_slack)
{
    constexpr int nx = static_cast<int>(NX);
    for(int k = 0; k <= Ni; ++k)
    {
        int xk_off = k * nx;
        for(int i = 0; i < nx; ++i)
            trips.emplace_back(row + i, xk_off + i, Scalar{1});
        if(need_slack)
        {
            int sk_off = n_states + k * nx;
            for(int i = 0; i < nx; ++i)
                trips.emplace_back(row + i, sk_off + i, Scalar{-1});
        }
        row += nx;
    }
}

/// Build MHE residual bound constraint triplets: H_k x_k
template <typename Scalar, std::size_t NX, std::size_t NY>
inline void add_mhe_residual_triplets(std::vector<Eigen::Triplet<Scalar>>& trips, int& row, int Ni, const Matrix<Scalar, NY, NX>& H)
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int ny = static_cast<int>(NY);
    for(int k = 0; k <= Ni; ++k)
    {
        int xk_off = k * nx;
        for(int i = 0; i < ny; ++i)
            for(int j = 0; j < nx; ++j)
                if(H(i, j) != Scalar{0})
                    trips.emplace_back(row + i, xk_off + j, H(i, j));
        row += ny;
    }
}

/// Build the MHE constraint matrix from dynamics, box, and residual blocks.
template <typename Scalar, std::size_t NX, std::size_t NY>
[[nodiscard]] auto build_mhe_constraint_matrix(const mhe_qp_dims& dims,
                                               int Ni,
                                               const Matrix<Scalar, NX, NX>& A,
                                               const Matrix<Scalar, NY, NX>& H,
                                               bool has_box_bounds,
                                               bool has_soft_constraints,
                                               bool has_residual_bounds) -> Eigen::SparseMatrix<Scalar, Eigen::ColMajor>
{
    constexpr int nx = static_cast<int>(NX);
    std::vector<Eigen::Triplet<Scalar>> a_trips;
    a_trips.reserve(static_cast<std::size_t>(Ni * (nx * nx + nx) + dims.n_box * 2 + dims.n_residual * nx));

    int row = 0;
    add_mhe_dynamics_triplets<Scalar, NX>(a_trips, row, Ni, A);
    if(has_box_bounds)
        add_mhe_box_triplets<Scalar, NX>(a_trips, row, Ni, dims.n_states, has_soft_constraints);
    if(has_residual_bounds)
        add_mhe_residual_triplets<Scalar, NX, NY>(a_trips, row, Ni, H);

    Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A_con(dims.n_con, dims.n_dec);
    A_con.setFromTriplets(a_trips.begin(), a_trips.end());
    return A_con;
}

/// Build the constant parts of the MHE QP: Hessian P and constraint matrix A.
/// @cite rao2003 -- Rao et al., "Constrained State Estimation for Nonlinear Discrete-Time Systems", 2003
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
[[nodiscard]] auto build_mhe_qp_structure(std::size_t N,
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
    int Ni = static_cast<int>(N);
    auto dims = compute_mhe_dims<NX, NY>(N, has_box_bounds, has_soft_constraints, has_residual_bounds);

    auto P = build_mhe_hessian<Scalar, NX, NY>(dims, Ni, arrival_weight, P_arr_inv, Q_inv, R_inv, A_lin[0], H_lin[0]);
    auto A_con = build_mhe_constraint_matrix<Scalar, NX, NY>(dims, Ni, A_lin[0], H_lin[0], has_box_bounds, has_soft_constraints, has_residual_bounds);

    Eigen::VectorX<Scalar> q = Eigen::VectorX<Scalar>::Zero(dims.n_dec);
    if(dims.n_slack > 0)
    {
        int slack_off = dims.n_states;
        for(int i = 0; i < dims.n_slack; ++i)
            q(slack_off + i) = soft_penalty;
    }

    Eigen::VectorX<Scalar> l = Eigen::VectorX<Scalar>::Zero(dims.n_con);
    Eigen::VectorX<Scalar> u = Eigen::VectorX<Scalar>::Zero(dims.n_con);

    return {std::move(P), std::move(q), std::move(A_con), std::move(l), std::move(u)};
}

/// Build the linear cost vector q for the MHE QP update.
/// @cite rao2003 -- Linear cost from arrival cost, process noise cross-terms, and measurement terms
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
[[nodiscard]] auto build_mhe_linear_cost(const mhe_qp_dims& dims,
                                         int Ni,
                                         Scalar arrival_weight,
                                         const Matrix<Scalar, NX, NX>& P_arr_inv,
                                         const Matrix<Scalar, NX, NX>& Q_inv,
                                         const Matrix<Scalar, NX, NX>& A_lin,
                                         const Matrix<Scalar, NX, NU>& B_lin,
                                         const Matrix<Scalar, NY, NX>& H_lin,
                                         const Matrix<Scalar, NY, NY>& R_inv,
                                         const Vector<Scalar, NX>& x_arrival,
                                         std::span<const Vector<Scalar, NU>> u_buf,
                                         std::span<const Vector<Scalar, NY>> z_buf,
                                         Scalar soft_penalty) -> Eigen::VectorX<Scalar>
{
    constexpr int nx = static_cast<int>(NX);

    Matrix<Scalar, NX, NX> AtQinv = (A_lin.transpose() * Q_inv).eval();
    Matrix<Scalar, NX, NY> HtRinv = (H_lin.transpose() * R_inv).eval();

    Eigen::VectorX<Scalar> q = Eigen::VectorX<Scalar>::Zero(dims.n_dec);
    q.template segment<nx>(0) = -arrival_weight * P_arr_inv * x_arrival;

    for(int k = 0; k < Ni; ++k)
    {
        Vector<Scalar, NX> Buk = B_lin * u_buf[static_cast<std::size_t>(k)];
        q.segment(k * nx, nx) += AtQinv * Buk;
        q.segment((k + 1) * nx, nx) -= Q_inv * Buk;
    }

    for(int k = 0; k <= Ni; ++k)
        q.segment(k * nx, nx) -= HtRinv * z_buf[static_cast<std::size_t>(k)];

    if(dims.n_slack > 0)
        for(int i = 0; i < dims.n_slack; ++i)
            q(dims.n_states + i) = soft_penalty;

    return q;
}

/// Build the MHE QP update bounds (l, u) from window data.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
[[nodiscard]] auto build_mhe_update_bounds(const mhe_qp_dims& dims,
                                           int Ni,
                                           const Matrix<Scalar, NX, NU>& B_lin,
                                           std::span<const Vector<Scalar, NU>> u_buf,
                                           std::span<const Vector<Scalar, NY>> z_buf,
                                           bool has_box_bounds,
                                           const std::optional<Vector<Scalar, NX>>& x_min,
                                           const std::optional<Vector<Scalar, NX>>& x_max,
                                           bool has_residual_bounds,
                                           const std::optional<Vector<Scalar, NY>>& residual_bound) -> std::pair<Eigen::VectorX<Scalar>, Eigen::VectorX<Scalar>>
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int ny = static_cast<int>(NY);
    constexpr auto inf = std::numeric_limits<Scalar>::infinity();

    Eigen::VectorX<Scalar> l(dims.n_con);
    Eigen::VectorX<Scalar> u(dims.n_con);
    int row = 0;

    for(int k = 0; k < Ni; ++k)
    {
        Vector<Scalar, NX> rhs = B_lin * u_buf[static_cast<std::size_t>(k)];
        l.segment(row, nx) = rhs;
        u.segment(row, nx) = rhs;
        row += nx;
    }

    if(has_box_bounds)
    {
        Vector<Scalar, NX> lb = x_min.value_or(Vector<Scalar, NX>::Constant(-inf));
        Vector<Scalar, NX> ub = x_max.value_or(Vector<Scalar, NX>::Constant(inf));
        for(int k = 0; k <= Ni; ++k)
        {
            l.segment(row, nx) = lb;
            u.segment(row, nx) = ub;
            row += nx;
        }
    }

    if(has_residual_bounds && residual_bound.has_value())
    {
        const auto& thresh = residual_bound.value();
        for(int k = 0; k <= Ni; ++k)
        {
            const auto& zk = z_buf[static_cast<std::size_t>(k)];
            for(int i = 0; i < ny; ++i)
            {
                l(row + i) = zk(i) - thresh(i);
                u(row + i) = zk(i) + thresh(i);
            }
            row += ny;
        }
    }

    return {std::move(l), std::move(u)};
}

/// Build per-solve QP update vectors (q, l, u) from current window data.
/// @cite rao2003 -- Called each time the MHE window shifts and a new solve is needed
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
[[nodiscard]] auto build_mhe_qp_update(std::size_t N,
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
    int Ni = static_cast<int>(N);
    auto dims = compute_mhe_dims<NX, NY>(N, has_box_bounds, has_soft_constraints, has_residual_bounds);

    auto q = build_mhe_linear_cost<Scalar, NX, NU, NY>(dims, Ni, arrival_weight, P_arr_inv, Q_inv, A_lin, B_lin, H_lin, R_inv, x_arrival, u_buf, z_buf, soft_penalty);
    auto [l, u] = build_mhe_update_bounds<Scalar, NX, NU, NY>(dims, Ni, B_lin, u_buf, z_buf, has_box_bounds, x_min, x_max, has_residual_bounds, residual_bound);

    return {std::move(q), std::move(l), std::move(u), warm_x, warm_y};
}

} // namespace ctrlpp::detail

#endif
