#ifndef HPP_GUARD_CTRLPP_CONTROL_DARE_H
#define HPP_GUARD_CTRLPP_CONTROL_DARE_H

#include "ctrlpp/types.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <complex>
#include <cstddef>
#include <optional>

namespace ctrlpp {

namespace detail {

// Swap two adjacent 1x1 blocks at positions p and p+1 in the upper triangular
// complex Schur form T, applying the corresponding unitary transformation to U.
//
// Given M = [[a, c], [0, b]], find unitary G such that G^H M G = [[b, c'], [0, a]].
// Method: solve a*x - x*b = -c => x = c/(b-a), then construct Givens from [x, 1].
template <typename Scalar, int N>
void swap_complex_schur_1x1(
    Eigen::Matrix<std::complex<Scalar>, N, N> &T,
    Eigen::Matrix<std::complex<Scalar>, N, N> &U,
    int p)
{
    using Complex = std::complex<Scalar>;

    Complex a = T(p, p);
    Complex b = T(p + 1, p + 1);

    if(std::abs(b - a) < Scalar{1e-14})
        return;

    Complex c = T(p, p + 1);

    Complex diff = b - a;
    Complex cs, sn;
    if(std::abs(c) < Scalar{1e-14})
    {
        // Off-diagonal is zero: use pure permutation (swap)
        cs = Complex{0};
        sn = Complex{1};
    }
    else
    {
        // Givens rotation that zeros (1,0) element after swap:
        // cs/sn = c/(b-a), so construct from vector [c, b-a]
        Scalar norm_val = std::sqrt(std::norm(c) + std::norm(diff));
        cs = c / norm_val;
        sn = diff / norm_val;
    }

    // G = [[cs, -conj(sn)], [sn, conj(cs)]]
    // T_new = G^H T G, U_new = U G

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

// Reorder a complex Schur decomposition so that stable eigenvalues (|lambda| < 1)
// appear in the top-left block. Returns the number of stable eigenvalues placed.
template <typename Scalar, int N>
auto reorder_complex_schur_stable_first(
    Eigen::Matrix<std::complex<Scalar>, N, N> &T,
    Eigen::Matrix<std::complex<Scalar>, N, N> &U,
    int required_stable) -> int
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

}

// Discrete Algebraic Riccati Equation solver using symplectic Schur decomposition.
// Solves: A^T P A - P - A^T P B (R + B^T P B)^{-1} B^T P A + Q = 0
// Returns the stabilizing solution P, or std::nullopt if the system is not stabilizable.
template <typename Scalar, std::size_t NX, std::size_t NU>
auto dare(const Eigen::Matrix<Scalar, int(NX), int(NX)> &A,
          const Eigen::Matrix<Scalar, int(NX), int(NU)> &B,
          const Eigen::Matrix<Scalar, int(NX), int(NX)> &Q,
          const Eigen::Matrix<Scalar, int(NU), int(NU)> &R)
    -> std::optional<Eigen::Matrix<Scalar, int(NX), int(NX)>>
{
    constexpr int n = static_cast<int>(NX);
    constexpr int n2 = 2 * n;

    using Complex = std::complex<Scalar>;
    using MatNxN = Eigen::Matrix<Scalar, n, n>;
    using Mat2Nx2N = Eigen::Matrix<Scalar, n2, n2>;
    using MatC2Nx2N = Eigen::Matrix<Complex, n2, n2>;
    using MatCNxN = Eigen::Matrix<Complex, n, n>;

    // Check A is invertible
    auto qr_At = A.transpose().colPivHouseholderQr();
    if(!qr_At.isInvertible())
        return std::nullopt;

    MatNxN I = MatNxN::Identity();
    MatNxN AinvT = qr_At.solve(I).eval();

    // G = B R^{-1} B^T
    MatNxN G = (B * R.colPivHouseholderQr().solve(
        Eigen::Matrix<Scalar, int(NU), int(NX)>(B.transpose())
    )).eval();

    // Form symplectic matrix Z (2n x 2n)
    Mat2Nx2N Z;
    Z.template block<n, n>(0, 0) = A + G * AinvT * Q;
    Z.template block<n, n>(0, n) = -G * AinvT;
    Z.template block<n, n>(n, 0) = -AinvT * Q;
    Z.template block<n, n>(n, n) = AinvT;

    // Complex Schur decomposition: Z = U T U^H
    Eigen::ComplexSchur<Mat2Nx2N> schur(Z);
    MatC2Nx2N T = schur.matrixT();
    MatC2Nx2N U = schur.matrixU();

    // Reorder so that stable eigenvalues (|lambda| < 1) are in top-left
    int stable = detail::reorder_complex_schur_stable_first(T, U, n);
    if(stable < n)
        return std::nullopt;

    // Extract U11 (top-left n x n) and U21 (bottom-left n x n) from reordered U
    MatCNxN U11 = U.template block<n, n>(0, 0);
    MatCNxN U21 = U.template block<n, n>(n, 0);

    // Check U11 is invertible
    auto qr_U11 = U11.colPivHouseholderQr();
    if(!qr_U11.isInvertible())
        return std::nullopt;

    // P = U21 * U11^{-1}
    MatCNxN Pc = (U21 * qr_U11.inverse()).eval();

    // Take real part (P should be real for a real-valued system)
    MatNxN P = Pc.real();

    // Symmetrize
    P = (Scalar{0.5} * (P + P.transpose())).eval();

    // Validate positive semi-definite
    Eigen::SelfAdjointEigenSolver<MatNxN> eigsolver(P, Eigen::EigenvaluesOnly);
    for(int i = 0; i < n; ++i)
    {
        if(eigsolver.eigenvalues()(i) < Scalar{-1e-10})
            return std::nullopt;
    }

    return P;
}

// DARE with cross-weight N:
// Transforms to standard form via Q' = Q - N R^{-1} N^T, A' = A - B R^{-1} N^T
template <typename Scalar, std::size_t NX, std::size_t NU>
std::optional<Eigen::Matrix<Scalar, int(NX), int(NX)>> dare(const Eigen::Matrix<Scalar, int(NX), int(NX)> &A, const Eigen::Matrix<Scalar, int(NX), int(NU)> &B,
                                                            const Eigen::Matrix<Scalar, int(NX), int(NX)> &Q, const Eigen::Matrix<Scalar, int(NU), int(NU)> &R,
                                                            const Eigen::Matrix<Scalar, int(NX), int(NU)> &N)
{
    auto Rinv_Nt = R.colPivHouseholderQr().solve(Eigen::Matrix<Scalar, int(NU), int(NX)>(N.transpose())).eval();

    Eigen::Matrix<Scalar, int(NX), int(NX)> Qp = (Q - N * Rinv_Nt).eval();
    Eigen::Matrix<Scalar, int(NX), int(NX)> Ap = (A - B * Rinv_Nt).eval();

    return dare<Scalar, NX, NU>(Ap, B, Qp, R);
}

}

#endif
