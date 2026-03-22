#ifndef HPP_GUARD_CTRLPP_SYSID_N4SID_H
#define HPP_GUARD_CTRLPP_SYSID_N4SID_H

#include "ctrlpp/sysid/fit_metrics.h"
#include "ctrlpp/sysid/sysid_result.h"
#include "ctrlpp/state_space.h"
#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cstddef>

namespace ctrlpp {

namespace detail {

template<typename Derived>
[[nodiscard]] auto block_hankel(const Eigen::MatrixBase<Derived>& signal,
                                Eigen::Index start,
                                Eigen::Index block_rows,
                                Eigen::Index cols) -> Eigen::MatrixX<typename Derived::Scalar>
{
    auto ny = signal.rows();
    Eigen::MatrixX<typename Derived::Scalar> H(ny * block_rows, cols);

    for (Eigen::Index bi = 0; bi < block_rows; ++bi) {
        for (Eigen::Index c = 0; c < cols; ++c) {
            H.block(bi * ny, c, ny, 1) = signal.col(start + bi + c);
        }
    }
    return H;
}

// N4SID core: build Hankel matrices and compute weighted oblique projection via LQ.
// Returns the matrix whose SVD reveals observability structure.
template<typename Scalar>
struct n4sid_lq_result {
    Eigen::MatrixX<Scalar> O_i;
    Eigen::MatrixX<Scalar> Y_future;
    Eigen::MatrixX<Scalar> U_future;
    Eigen::Index i;
    Eigen::Index j;
};

template<typename Derived1, typename Derived2>
[[nodiscard]] auto n4sid_oblique_projection(const Eigen::MatrixBase<Derived1>& Y,
                                             const Eigen::MatrixBase<Derived2>& U,
                                             Eigen::Index i)
    -> n4sid_lq_result<typename Derived1::Scalar>
{
    using Scalar = typename Derived1::Scalar;

    auto N = Y.cols();
    auto ny = Y.rows();
    auto nu = U.rows();
    Eigen::Index j = N - 2 * i + 1;

    // Build block Hankel matrices
    auto Y_past   = block_hankel(Y, 0, i, j);
    auto Y_future = block_hankel(Y, i, i, j);
    auto U_past   = block_hankel(U, 0, i, j);
    auto U_future = block_hankel(U, i, i, j);

    // Stack: H = [U_f; W_p; Y_f] where W_p = [U_past; Y_past]
    // Row partition sizes:
    auto r1 = nu * i;   // U_future rows
    auto r2 = (nu + ny) * i;  // W_p rows
    auto r3 = ny * i;   // Y_future rows

    Eigen::MatrixX<Scalar> H(r1 + r2 + r3, j);
    H.topRows(r1) = U_future;
    H.middleRows(r1, nu * i) = U_past;
    H.middleRows(r1 + nu * i, ny * i) = Y_past;
    H.bottomRows(r3) = Y_future;

    // LQ decomposition: H = L * Q where L is lower triangular
    // Eigen provides HouseholderQR on H^T: H^T = Q_h * R_h, so H = R_h^T * Q_h^T = L * Q
    Eigen::HouseholderQR<Eigen::MatrixX<Scalar>> qr(H.transpose());
    auto total_rows = r1 + r2 + r3;
    auto r_size = std::min(j, total_rows);
    Eigen::MatrixX<Scalar> R_h = qr.matrixQR().topRows(r_size).template triangularView<Eigen::Upper>();
    Eigen::MatrixX<Scalar> L = R_h.transpose();  // L is lower triangular (total_rows x r_size)

    // Extract L blocks according to row partitioning:
    // L11: rows 0..r1-1, cols 0..r1-1           (U_f on U_f)
    // L21: rows r1..r1+r2-1, cols 0..r1-1       (W_p on U_f)
    // L22: rows r1..r1+r2-1, cols r1..r1+r2-1   (W_p on W_p)
    // L31: rows r1+r2..end, cols 0..r1-1         (Y_f on U_f)
    // L32: rows r1+r2..end, cols r1..r1+r2-1     (Y_f on W_p)
    // L33: rows r1+r2..end, cols r1+r2..end      (Y_f residual)

    // The oblique projection O_i = Y_f / U_f onto W_p is:
    // O_i = L32 * Q2 (where Q2 spans the W_p-related columns of Q)
    // For the SVD, we only need L32 since Q2 is orthonormal and doesn't affect singular values.

    auto l_rows = L.rows();
    auto l_cols = L.cols();

    // Clamp block sizes to available dimensions
    auto c1 = std::min(r1, l_cols);
    auto c2 = std::min(r2, l_cols - c1);
    auto br3 = std::min(r3, l_rows - r1 - r2);

    if (c2 <= 0 || br3 <= 0) {
        return {Eigen::MatrixX<Scalar>::Zero(r3, 1), Y_future, U_future, i, j};
    }

    Eigen::MatrixX<Scalar> L32 = L.block(r1 + r2, c1, br3, c2);

    return {L32, Y_future, U_future, i, j};
}

}

template<typename Derived1, typename Derived2>
[[nodiscard]] auto n4sid_singular_values(const Eigen::MatrixBase<Derived1>& Y,
                                         const Eigen::MatrixBase<Derived2>& U,
                                         std::size_t block_rows = 0)
    -> Eigen::VectorX<typename Derived1::Scalar>
{
    using Scalar = typename Derived1::Scalar;

    auto N = Y.cols();
    auto i = static_cast<Eigen::Index>(block_rows);
    if (i == 0) {
        i = std::min(static_cast<Eigen::Index>(N / 4), Eigen::Index{30});
    }

    Eigen::Index j = N - 2 * i + 1;
    if (j <= 0) {
        return Eigen::VectorX<Scalar>{};
    }

    auto [O_i, Y_f, U_f, bi, bj] = detail::n4sid_oblique_projection(Y, U, i);

    Eigen::BDCSVD<Eigen::MatrixX<Scalar>> svd(O_i, Eigen::ComputeThinU | Eigen::ComputeThinV);

    return svd.singularValues();
}

template<std::size_t NX, typename Derived1, typename Derived2>
[[nodiscard]] auto n4sid(const Eigen::MatrixBase<Derived1>& Y,
                         const Eigen::MatrixBase<Derived2>& U,
                         std::size_t block_rows = 0)
    -> n4sid_result<typename Derived1::Scalar, NX, 1, 1>
{
    using Scalar = typename Derived1::Scalar;
    static constexpr auto nx = static_cast<Eigen::Index>(NX);

    auto N = Y.cols();
    auto ny = Y.rows();
    auto nu = U.rows();

    auto i = static_cast<Eigen::Index>(block_rows);
    if (i == 0) {
        i = std::min(static_cast<Eigen::Index>(N / 4), Eigen::Index{30});
    }

    auto [O_i, Y_future, U_future, bi, bj] = detail::n4sid_oblique_projection(Y, U, i);

    // BDCSVD of the oblique projection
    Eigen::BDCSVD<Eigen::MatrixX<Scalar>> svd(O_i, Eigen::ComputeThinU | Eigen::ComputeThinV);

    auto sv = svd.singularValues();
    auto U_svd = svd.matrixU();

    // Condition number from the NX singular values
    Scalar cond = (sv.size() >= nx && sv(nx - 1) > Scalar{0})
        ? sv(0) / sv(nx - 1)
        : std::numeric_limits<Scalar>::infinity();

    // Observability matrix Gamma: first NX columns of U_svd * diag(sqrt(S))
    Eigen::Index rank = std::min(static_cast<Eigen::Index>(sv.size()), nx);
    Eigen::MatrixX<Scalar> Gamma(U_svd.rows(), nx);
    Gamma.setZero();
    for (Eigen::Index c = 0; c < rank; ++c) {
        Gamma.col(c) = U_svd.col(c) * std::sqrt(sv(c));
    }

    // C from first ny rows of Gamma
    Matrix<Scalar, 1, NX> C;
    for (Eigen::Index c = 0; c < nx; ++c) {
        C(0, c) = Gamma(c);
    }

    // A from shift relation on observability matrix:
    // Gamma = [C; CA; CA^2; ...; CA^{i-1}]
    // Gamma_shifted = Gamma without first ny rows = [CA; CA^2; ...; CA^{i-1}]
    // Gamma_trunc   = Gamma without last ny rows  = [C; CA; ...; CA^{i-2}]
    // Relation: Gamma_shifted = Gamma_trunc * A
    // Solve: A = Gamma_trunc \ Gamma_shifted
    auto gamma_rows = Gamma.rows();
    Eigen::MatrixX<Scalar> Gamma_trunc = Gamma.topRows(gamma_rows - ny);
    Eigen::MatrixX<Scalar> Gamma_shifted = Gamma.bottomRows(gamma_rows - ny);

    Eigen::MatrixX<Scalar> A_dyn = Gamma_trunc.colPivHouseholderQr().solve(Gamma_shifted);

    Matrix<Scalar, NX, NX> A;
    for (Eigen::Index r = 0; r < nx; ++r) {
        for (Eigen::Index c = 0; c < nx; ++c) {
            A(r, c) = A_dyn(r, c);
        }
    }

    // Recover B and D from original time-domain data via least squares.
    // With A, C known, solve for B, D that minimize simulation error.
    // For each time step: y(t) = C*x(t) + D*u(t), x(t+1) = A*x(t) + B*u(t)
    //
    // Joint least-squares: stack the system equations.
    // For SISO: B is (NX x 1), D is (1 x 1), so we have NX+1 unknowns.
    // From the state equation:  x(t+1) = A*x(t) + B*u(t)
    // From the output equation: y(t) = C*x(t) + D*u(t)
    //
    // Build a regressor: for each t, the contribution of [B; D] to
    // [x(t+1); y(t)] given x(t) and u(t).
    //
    // Efficient approach: simulate with A,C to get state sequence given
    // candidate B, then solve the linear system for B, D.

    // Method: iterative refinement starting from zero B, D
    // Actually, direct approach: write out the dependency and solve via
    // the input-to-state transfer relation.
    //
    // x(t) = A^t * x(0) + sum_{k=0}^{t-1} A^{t-1-k} * B * u(k)
    // y(t) = C * x(t) + D * u(t)
    //
    // With x(0) = 0:
    // y(t) = sum_{k=0}^{t-1} C * A^{t-1-k} * B * u(k) + D * u(t)
    //       = sum_{k=0}^{t-1} (C * A^{t-1-k}) * B * u(k) + D * u(t)
    //
    // This is linear in [B; D] (NX+1 unknowns for SISO).
    // Build regressor Phi: each row corresponds to one time step.
    // Phi(t, :) = [sum_{k=0}^{t-1} C*A^{t-1-k}*u(k) (vectorized as NX entries), u(t)]

    // For SISO: B is (NX x 1), so Phi_B(t) is a (1 x NX) vector
    // Phi_B(t, j) = sum_{k=0}^{t-1} [C*A^{t-1-k}]_j * u(k) for the j-th element of B
    // Phi_D(t) = u(t)

    // Actually simpler: precompute Markov parameters h(k) = C*A^k for k=0..N-1
    // Then y(t) = sum_{k=0}^{t-1} h(t-1-k) * B * u(k) + D*u(t)
    // = [conv(h, u)](t) * B + D*u(t)   (but h is a (1 x NX) row, B is (NX x 1))

    // Build Phi matrix for least squares
    // Unknowns: bd = [B(0); B(1); ...; B(NX-1); D]   (NX+1 x 1)
    Eigen::MatrixX<Scalar> Phi(N, nx + nu);
    Phi.setZero();
    Eigen::VectorX<Scalar> y_target(N);

    // Precompute powers of A applied through C
    // For each time t, compute the convolution contribution
    // x_contrib(j, t) = sum_{k=0}^{t-1} [A^{t-1-k}](:,j) * u(k)
    // Actually: maintain x_basis(j) = e_j, propagate each through A
    // More efficient: maintain NX "basis state" vectors, each propagated by A

    // State basis propagation: x_j(t+1) = A * x_j(t) + e_j * u(t)
    // Then Phi_B(t, j) = C * x_j(t) and Phi_D(t) = u(t)
    Eigen::MatrixX<Scalar> x_basis = Eigen::MatrixX<Scalar>::Zero(nx, nx);  // columns = basis states

    for (Eigen::Index t = 0; t < N; ++t) {
        Scalar u_t = U(0, t);
        y_target(t) = Y(0, t);

        // Phi row: [C * x_basis(:, 0), C * x_basis(:, 1), ..., u(t)]
        for (Eigen::Index j = 0; j < nx; ++j) {
            Phi(t, j) = (C * x_basis.col(j))(0, 0);
        }
        Phi(t, nx) = u_t;

        // Propagate all basis states: x_j(t+1) = A * x_j(t) + e_j * u(t)
        Eigen::MatrixX<Scalar> x_basis_new = A_dyn * x_basis;
        for (Eigen::Index j = 0; j < nx; ++j) {
            x_basis_new(j, j) += u_t;
        }
        x_basis = x_basis_new;
    }

    // Solve: y_target = Phi * bd
    Eigen::VectorX<Scalar> bd = Phi.colPivHouseholderQr().solve(y_target);

    Matrix<Scalar, NX, 1> B;
    for (Eigen::Index r = 0; r < nx; ++r) {
        B(r, 0) = bd(r);
    }

    Matrix<Scalar, 1, 1> D;
    D(0, 0) = bd(nx);

    discrete_state_space<Scalar, NX, 1, 1> sys{.A = A, .B = B, .C = C, .D = D};

    // Simulate identified model on original data for fit metrics
    Vector<Scalar, NX> x = Vector<Scalar, NX>::Zero();
    Eigen::VectorX<Scalar> y_predicted(N);
    Eigen::VectorX<Scalar> y_actual(N);

    for (Eigen::Index t = 0; t < N; ++t) {
        Matrix<Scalar, 1, 1> u_vec;
        u_vec(0, 0) = U(0, t);
        auto y_hat = (C * x + D * u_vec).eval();
        y_predicted(t) = y_hat(0, 0);
        x = (A * x + B * u_vec).eval();
        y_actual(t) = Y(0, t);
    }

    auto metrics = compute_fit_metrics(y_actual, y_predicted);

    return {.system = sys, .singular_values = sv, .metrics = metrics, .condition_number = cond};
}

}

#endif
