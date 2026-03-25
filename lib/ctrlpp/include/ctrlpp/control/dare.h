#ifndef HPP_GUARD_CTRLPP_CONTROL_DARE_H
#define HPP_GUARD_CTRLPP_CONTROL_DARE_H

/// @brief Discrete Algebraic Riccati Equation solver via symplectic Schur decomposition.
///
/// @cite laub1979 -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979

#include "ctrlpp/types.h"

#include "ctrlpp/detail/covariance_ops.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <complex>
#include <cstddef>
#include <optional>

namespace ctrlpp
{

namespace detail
{

// Compute Givens rotation parameters for swapping adjacent 1x1 Schur blocks.
template <typename Scalar>
auto compute_schur_swap_givens(std::complex<Scalar> a, std::complex<Scalar> b, std::complex<Scalar> c)
    -> std::pair<std::complex<Scalar>, std::complex<Scalar>>
{
    using Complex = std::complex<Scalar>;
    Complex diff = b - a;
    if(std::abs(c) < Scalar{1e-14})
        return {Complex{0}, Complex{1}};

    Scalar norm_val = std::sqrt(std::norm(c) + std::norm(diff));
    return {c / norm_val, diff / norm_val};
}

// Apply Givens rotation G = [[cs, -conj(sn)], [sn, conj(cs)]] to T and U at position p.
template <typename Scalar, int N>
void apply_schur_givens(Eigen::Matrix<std::complex<Scalar>, N, N>& T,
                        Eigen::Matrix<std::complex<Scalar>, N, N>& U,
                        int p, std::complex<Scalar> cs, std::complex<Scalar> sn)
{
    using Complex = std::complex<Scalar>;

    // Left multiply rows p, p+1 by G^H:
    for(int j = 0; j < N; ++j)
    {
        Complex t0 = T(p, j);
        Complex t1 = T(p + 1, j);
        T(p, j) = std::conj(cs) * t0 + std::conj(sn) * t1;
        T(p + 1, j) = -sn * t0 + cs * t1;
    }
    // Right multiply cols p, p+1 by G:
    for(int i = 0; i < N; ++i)
    {
        Complex t0 = T(i, p);
        Complex t1 = T(i, p + 1);
        T(i, p) = t0 * cs + t1 * sn;
        T(i, p + 1) = t0 * (-std::conj(sn)) + t1 * std::conj(cs);
    }
    T(p + 1, p) = Complex{0};

    // Right multiply U by G:
    for(int i = 0; i < N; ++i)
    {
        Complex u0 = U(i, p);
        Complex u1 = U(i, p + 1);
        U(i, p) = u0 * cs + u1 * sn;
        U(i, p + 1) = u0 * (-std::conj(sn)) + u1 * std::conj(cs);
    }
}

// Swap two adjacent 1x1 blocks at positions p and p+1 in the upper triangular
// complex Schur form T, applying the corresponding unitary transformation to U.
/// @cite laub1979 -- Laub, 1979 (Schur reordering for DARE)
template <typename Scalar, int N>
void swap_complex_schur_1x1(Eigen::Matrix<std::complex<Scalar>, N, N>& T, Eigen::Matrix<std::complex<Scalar>, N, N>& U, int p)
{
    if(std::abs(T(p + 1, p + 1) - T(p, p)) < Scalar{1e-14})
        return;

    auto [cs, sn] = compute_schur_swap_givens(T(p, p), T(p + 1, p + 1), T(p, p + 1));
    apply_schur_givens(T, U, p, cs, sn);
}

// Reorder a complex Schur decomposition so that stable eigenvalues (|lambda| < 1)
// appear in the top-left block. Returns the number of stable eigenvalues placed.
template <typename Scalar, int N>
auto reorder_complex_schur_stable_first(Eigen::Matrix<std::complex<Scalar>, N, N>& T, Eigen::Matrix<std::complex<Scalar>, N, N>& U, int required_stable) -> int
{
    int stable_count = 0;

    while(stable_count < required_stable)
    {
        int pos = stable_count;
        bool found = false;

        while(pos < N)
        {
            if(std::abs(T(pos, pos)) < Scalar{1})
            {
                found = true;
                break;
            }
            ++pos;
        }

        if(!found)
            return stable_count;

        // Bubble from pos to stable_count via adjacent swaps
        for(int k = pos; k > stable_count; --k)
            swap_complex_schur_1x1(T, U, k - 1);

        ++stable_count;
    }

    return stable_count;
}

} // namespace detail

/// @brief Build symplectic matrix Z for DARE from system matrices A, B, Q, R.
///
/// @cite laub1979 -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979, Eq. 7
template <typename Scalar, std::size_t NX, std::size_t NU>
auto build_dare_symplectic(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
                           const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
                           const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
                           const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
    -> std::optional<Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>>
{
    constexpr int n = static_cast<int>(NX);
    constexpr int n2 = 2 * n;
    using MatNxN = Eigen::Matrix<Scalar, n, n>;
    using Mat2Nx2N = Eigen::Matrix<Scalar, n2, n2>;

    auto qr_At = A.transpose().colPivHouseholderQr();
    if(!qr_At.isInvertible())
        return std::nullopt;

    MatNxN AinvT = qr_At.solve(MatNxN::Identity()).eval();
    MatNxN G = (B * R.colPivHouseholderQr().solve(Eigen::Matrix<Scalar, int(NU), int(NX)>(B.transpose()))).eval();

    Mat2Nx2N Z;
    Z.template block<n, n>(0, 0) = A + G * AinvT * Q;
    Z.template block<n, n>(0, n) = -G * AinvT;
    Z.template block<n, n>(n, 0) = -AinvT * Q;
    Z.template block<n, n>(n, n) = AinvT;
    return Z;
}

/// @brief Extract DARE solution P from reordered Schur decomposition: P = real(U21 * U11^{-1}).
///
/// @cite laub1979 -- Laub, 1979, Eq. 10
template <typename Scalar, int N2>
auto extract_dare_solution(const Eigen::Matrix<std::complex<Scalar>, N2, N2>& U)
    -> std::optional<Eigen::Matrix<Scalar, N2 / 2, N2 / 2>>
{
    constexpr int n = N2 / 2;
    using MatNxN = Eigen::Matrix<Scalar, n, n>;
    using MatCNxN = Eigen::Matrix<std::complex<Scalar>, n, n>;

    MatCNxN U11 = U.template block<n, n>(0, 0);
    MatCNxN U21 = U.template block<n, n>(n, 0);

    auto qr_U11 = U11.colPivHouseholderQr();
    if(!qr_U11.isInvertible())
        return std::nullopt;

    MatNxN P = detail::symmetrize((U21 * qr_U11.inverse()).eval().real());

    Eigen::SelfAdjointEigenSolver<MatNxN> eigsolver(P, Eigen::EigenvaluesOnly);
    for(int i = 0; i < n; ++i)
    {
        if(eigsolver.eigenvalues()(i) < Scalar{-1e-10})
            return std::nullopt;
    }
    return P;
}

// Discrete Algebraic Riccati Equation solver using symplectic Schur decomposition.
// Solves: A^T P A - P - A^T P B (R + B^T P B)^{-1} B^T P A + Q = 0
// Returns the stabilizing solution P, or std::nullopt if the system is not stabilizable.
template <typename Scalar, std::size_t NX, std::size_t NU>
auto dare(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A, const Eigen::Matrix<Scalar, int(NX), int(NU)>& B, const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q, const Eigen::Matrix<Scalar, int(NU), int(NU)>& R)
    -> std::optional<Eigen::Matrix<Scalar, int(NX), int(NX)>>
{
    constexpr int n = static_cast<int>(NX);
    constexpr int n2 = 2 * n;

    auto Z_opt = build_dare_symplectic<Scalar, NX, NU>(A, B, Q, R);
    if(!Z_opt)
        return std::nullopt;

    Eigen::ComplexSchur<Eigen::Matrix<Scalar, n2, n2>> schur(*Z_opt);
    Eigen::Matrix<std::complex<Scalar>, n2, n2> T = schur.matrixT();
    Eigen::Matrix<std::complex<Scalar>, n2, n2> U = schur.matrixU();

    int stable = detail::reorder_complex_schur_stable_first(T, U, n);
    if(stable < n)
        return std::nullopt;

    return extract_dare_solution<Scalar, n2>(U);
}

// DARE with cross-weight N:
// Transforms to standard form via Q' = Q - N R^{-1} N^T, A' = A - B R^{-1} N^T
template <typename Scalar, std::size_t NX, std::size_t NU>
std::optional<Eigen::Matrix<Scalar, int(NX), int(NX)>> dare(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
                                                            const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
                                                            const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
                                                            const Eigen::Matrix<Scalar, int(NU), int(NU)>& R,
                                                            const Eigen::Matrix<Scalar, int(NX), int(NU)>& N)
{
    auto Rinv_Nt = R.colPivHouseholderQr().solve(Eigen::Matrix<Scalar, int(NU), int(NX)>(N.transpose())).eval();

    Eigen::Matrix<Scalar, int(NX), int(NX)> Qp = (Q - N * Rinv_Nt).eval();
    Eigen::Matrix<Scalar, int(NX), int(NX)> Ap = (A - B * Rinv_Nt).eval();

    return dare<Scalar, NX, NU>(Ap, B, Qp, R);
}

} // namespace ctrlpp

#endif
