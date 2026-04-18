#ifndef HPP_GUARD_CTRLPP_CONTROL_CARE_H
#define HPP_GUARD_CTRLPP_CONTROL_CARE_H

/// @brief Continuous-time Algebraic Riccati Equation solver via Hamiltonian Schur decomposition.
///
/// Solves A^T P + P A - P B R^{-1} B^T P + Q = 0 for the stabilizing P.
///
/// Builds the 2n x 2n Hamiltonian H = [[A, -B R^{-1} B^T], [-Q, -A^T]], computes its
/// complex Schur decomposition, reorders eigenvalues with Re(lambda) < 0 (continuous-stable
/// modes) to the top-left via a hand-rolled Givens-rotation bubble sort, and extracts
/// P = real(U21 * U11^{-1}) from the corresponding invariant subspace basis.
///
/// Reuses the complex-Schur reordering primitives and the P-extraction helper from
/// dare.h -- the continuous-time variant only differs in the Hamiltonian build and the
/// stability criterion applied per eigenvalue.
///
/// @cite laub1979 -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979

#include "ctrlpp/control/dare.h"
#include "ctrlpp/types.h"

#include "ctrlpp/detail/covariance_ops.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <complex>
#include <cstddef>
#include <optional>
#include <type_traits>

namespace ctrlpp
{

namespace detail
{

/// @brief Build the 2n x 2n Hamiltonian matrix H for CARE from (A, B, Q, R).
///
/// H = [[A,    -B R^{-1} B^T],
///      [-Q,   -A^T         ]]
template <typename Scalar, std::size_t NX, std::size_t NU>
auto build_care_hamiltonian(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
                            const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
                            const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
                            const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
    -> std::optional<Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>>
{
    constexpr int n = static_cast<int>(NX);
    constexpr int n2 = 2 * n;
    using MatNxN = Eigen::Matrix<Scalar, n, n>;
    using Mat2Nx2N = Eigen::Matrix<Scalar, n2, n2>;

    MatNxN S = (B * R.colPivHouseholderQr().solve(
                        Eigen::Matrix<Scalar, int(NU), int(NX)>(B.transpose())))
                   .eval();

    Mat2Nx2N H;
    H.template block<n, n>(0, 0) = A;
    H.template block<n, n>(0, n) = -S;
    H.template block<n, n>(n, 0) = -Q;
    H.template block<n, n>(n, n) = -A.transpose();

    if(!H.allFinite())
        return std::nullopt;

    return H;
}

/// @brief Reorder a complex Schur decomposition so eigenvalues with Re(lambda) < 0
/// appear in the top-left block. Returns the number of stable eigenvalues placed.
///
/// Mirrors `reorder_complex_schur_stable_first` from dare.h but applies the continuous-time
/// stability criterion instead of the discrete-time one.
template <typename Scalar, int N>
auto reorder_complex_schur_lhp_first(Eigen::Matrix<std::complex<Scalar>, N, N>& T,
                                     Eigen::Matrix<std::complex<Scalar>, N, N>& U,
                                     int required_stable) -> int
{
    int stable_count = 0;

    while(stable_count < required_stable)
    {
        int pos = stable_count;
        bool found = false;

        while(pos < N)
        {
            if(T(pos, pos).real() < Scalar{0})
            {
                found = true;
                break;
            }
            ++pos;
        }

        if(!found)
            return stable_count;

        for(int k = pos; k > stable_count; --k)
            swap_complex_schur_1x1<Scalar, N>(T, U, k - 1);

        ++stable_count;
    }

    return stable_count;
}

}

/// @brief Continuous-time Algebraic Riccati Equation solver using Hamiltonian Schur decomposition.
///
/// Solves A^T P + P A - P B R^{-1} B^T P + Q = 0 for the stabilizing P.
/// Returns std::nullopt if the Hamiltonian build fails, fewer than n eigenvalues lie in
/// the open left half-plane, or if the U11 Schur block is singular.
///
/// The continuous-time LQR gain is K = R^{-1} B^T P.
template <typename Scalar, std::size_t NX, std::size_t NU>
auto care(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
          const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
          const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
          const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
    -> std::optional<Eigen::Matrix<Scalar, int(NX), int(NX)>>
{
    static_assert(std::is_floating_point_v<Scalar>, "Scalar must be a floating-point type");
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");

    constexpr int n = static_cast<int>(NX);
    constexpr int n2 = 2 * n;

    auto H_opt = detail::build_care_hamiltonian<Scalar, NX, NU>(A, B, Q, R);
    if(!H_opt)
        return std::nullopt;

    Eigen::ComplexSchur<Eigen::Matrix<Scalar, n2, n2>> schur(*H_opt);
    if(schur.info() != Eigen::Success)
        return std::nullopt;

    Eigen::Matrix<std::complex<Scalar>, n2, n2> T = schur.matrixT();
    Eigen::Matrix<std::complex<Scalar>, n2, n2> U = schur.matrixU();

    if(!T.allFinite() || !U.allFinite())
        return std::nullopt;

    int stable = detail::reorder_complex_schur_lhp_first<Scalar, n2>(T, U, n);
    if(stable < n)
        return std::nullopt;

    if(!U.allFinite())
        return std::nullopt;

    return extract_dare_solution<Scalar, n2>(U);
}

/// @brief CARE with cross-weight N. Reduces to standard form via
/// Q' = Q - N R^{-1} N^T, A' = A - B R^{-1} N^T, then forwards to the standard CARE solver.
template <typename Scalar, std::size_t NX, std::size_t NU>
auto care(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
          const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
          const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
          const Eigen::Matrix<Scalar, int(NU), int(NU)>& R,
          const Eigen::Matrix<Scalar, int(NX), int(NU)>& N)
    -> std::optional<Eigen::Matrix<Scalar, int(NX), int(NX)>>
{
    auto Rinv_Nt = R.colPivHouseholderQr()
                       .solve(Eigen::Matrix<Scalar, int(NU), int(NX)>(N.transpose()))
                       .eval();

    Eigen::Matrix<Scalar, int(NX), int(NX)> Qp = (Q - N * Rinv_Nt).eval();
    Eigen::Matrix<Scalar, int(NX), int(NX)> Ap = (A - B * Rinv_Nt).eval();

    return care<Scalar, NX, NU>(Ap, B, Qp, R);
}

}

#endif
