#ifndef HPP_GUARD_CTRLPP_DETAIL_COVARIANCE_OPS_H
#define HPP_GUARD_CTRLPP_DETAIL_COVARIANCE_OPS_H

/// @brief Covariance-matrix utilities: the counted rounding bound, admission of
/// a covariance (drift symmetrized, structural asymmetry and indefiniteness
/// refused), symmetrization, the matrix square root and the one-output gain.
///
/// @cite simon2006 -- Simon, "Optimal State Estimation", 2006, Ch. 5 (covariance conditioning / Joseph form)
/// @cite higham1988 -- Higham, "Computing a nearest symmetric positive semidefinite matrix", Linear Algebra Appl. 103, 1988
/// @cite higham2002 -- Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed., 2002, Sec. 2.1 and 3.1
/// @cite golub2013 -- Golub & Van Loan, "Matrix Computations", 4th ed., 2013, Sec. 8.1 (symmetric eigenvalue perturbation)

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

#include <limits>
#include <cstdint>
#include <algorithm>

namespace ctrlpp::detail
{

/// @brief The rounding bound `factor * count * eps * scale`, never below the
/// underflow term `count * denorm_min`.
///
/// `count` is in machine epsilons, eps = 2 u, derived by each caller from the
/// convention that an n-term inner product carries at most gamma_n, about n u,
/// of its scale (@cite higham2002, Sec. 3.1); `factor` is the order-one
/// prefactor left to the caller. The absolute term is the gradual-underflow
/// part of the standard model, fl(x op y) = (x op y)(1 + delta) + eta with
/// |eta| <= denorm_min / 2 (@cite higham2002, Sec. 2.1): without it the bound
/// collapses to zero below scale = denorm_min / (factor eps), and at and above
/// that scale the bound is the relative term exactly.
template <typename Scalar>
inline auto rounding_bound(std::int32_t count, Scalar scale, Scalar factor) -> Scalar
{
    const Scalar n = static_cast<Scalar>(count);
    return std::max(factor * n * std::numeric_limits<Scalar>::epsilon() * scale, n * std::numeric_limits<Scalar>::denorm_min());
}

/// @brief Symmetrize M in place when its asymmetry is rounding drift, and
/// refuse it, unchanged, when the asymmetry is structural.
///
/// The asymmetry max|M_ij - M_ji| is drift when it is within
/// `rounding_bound(N, max|M_ij|, factor)`. The count behind `factor = 1`: each
/// triangle of a matrix formed by products carries at most N u of its scale,
/// so the two triangles of one exactly symmetric matrix differ by at most
/// 2 N u = N eps times it. Each differing pair is replaced in both positions by
/// its midpoint M_ij + (M_ji - M_ij) / 2, which cannot overflow; an equal pair
/// is never written, so a symmetric M is left bitwise untouched. A non-finite
/// M is refused: an infinite entry would make the bound infinite too.
template <typename Scalar, int N>
inline auto symmetrize_within_drift(Eigen::Matrix<Scalar, N, N>& M, Scalar factor = Scalar{1}) -> bool
{
    if(!M.allFinite())
        return false;
    const Scalar bound = rounding_bound(N, M.cwiseAbs().maxCoeff(), factor);
    if(!((M - M.transpose()).cwiseAbs().array() <= bound).all())
        return false;
    for(int j = 1; j < N; ++j)
        for(int i = 0; i < j; ++i)
            if(M(i, j) != M(j, i))
                M(i, j) = M(j, i) = M(i, j) + (M(j, i) - M(i, j)) / Scalar{2};
    return true;
}

/// @brief Symmetrize a square matrix: 0.5 * (M + M^T)
///
/// Numerical operations on covariance matrices can introduce small asymmetries.
/// This enforces exact symmetry, which downstream Cholesky/LDLT decompositions require.
///
/// @cite simon2006 -- Simon, "Optimal State Estimation", 2006, Ch. 5 (numerical covariance conditioning)
template <typename Derived>
inline auto symmetrize(const Eigen::MatrixBase<Derived>& M)
{
    using Scalar = typename Derived::Scalar;
    return (Scalar{0.5} * (M + M.transpose())).eval();
}

/// @brief Result of a covariance square root: a factor S with S*S^T equal to
/// the (possibly repaired) covariance, and a flag recording whether the input
/// had to be repaired to the nearest symmetric positive definite matrix.
template <typename Scalar, int N>
struct covariance_sqrt_result
{
    Eigen::Matrix<Scalar, N, N> factor;
    bool repaired;
};

/// @brief Matrix square root S of a covariance P such that S*S^T = P.
///
/// The primary factor is the unpivoted lower Cholesky factor L (P = L*L^T),
/// which squares back to P exactly with no permutation to track. A pivoted
/// LDLT instead satisfies Pi*P*Pi^T = L*D*L^T for a permutation Pi, so reading
/// L*sqrt(D) off it reconstructs a symmetric permutation of P rather than P
/// itself; the unpivoted Cholesky avoids that trap entirely.
///
/// When P is not positive definite (Cholesky fails), P is repaired to the
/// nearest symmetric positive definite matrix in the Frobenius sense by a
/// symmetric eigendecomposition whose eigenvalues are clamped up to a
/// scale-aware floor, and the factor is rebuilt from the clamped spectrum. The
/// returned `repaired` flag records this so the caller can surface a degraded
/// filter-health status rather than continuing silently.
///
/// The clamp floor is `floor_factor * N * eps * max|P_ij|`. The `N * eps *
/// max|P_ij|` part is the order of the symmetric eigensolver's backward error
/// (Weyl perturbation bound), below which a computed eigenvalue is
/// indistinguishable from zero and must be lifted to keep the repaired matrix
/// strictly positive definite. `floor_factor` is the leading order-one
/// prefactor of that bound, exposed as a parameter so a caller may tighten or
/// loosen the floor for a particular problem.
///
/// @cite higham1988 -- Higham, "Computing a nearest symmetric positive semidefinite matrix", 1988
/// @cite golub2013 -- Golub & Van Loan, "Matrix Computations", 4th ed., 2013, Sec. 8.1
template <typename Scalar, int N>
inline covariance_sqrt_result<Scalar, N> covariance_sqrt(const Eigen::Matrix<Scalar, N, N>& P, Scalar floor_factor = Scalar{1})
{
    Eigen::LLT<Eigen::Matrix<Scalar, N, N>> llt(P);
    if(llt.info() == Eigen::Success)
        return {llt.matrixL(), false};

    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    const Scalar floor = floor_factor * static_cast<Scalar>(N) * eps * P.cwiseAbs().maxCoeff();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, N, N>> eig(symmetrize(P));
    Eigen::Matrix<Scalar, N, 1> clamped = eig.eigenvalues().cwiseMax(floor);
    Eigen::Matrix<Scalar, N, N> factor = eig.eigenvectors() * clamped.cwiseSqrt().asDiagonal();
    return {factor, true};
}

/// @brief Symmetric eigendecomposition of a candidate covariance and the
/// verdict on whether it is positive semidefinite.
template <typename Scalar, int N>
struct psd_decomposition_result
{
    Eigen::Matrix<Scalar, N, 1> eigenvalues;
    Eigen::Matrix<Scalar, N, N> eigenvectors;
    bool psd;
};

/// @brief Decide whether M is positive semidefinite, refusing rather than
/// repairing.
///
/// M is PSD when its smallest eigenvalue is at least `-floor`, with the floor
/// `rounding_bound(N, max|M_ij|, floor_factor)`. Its relative term is the
/// `floor_factor * N * eps * max|M_ij|` that `covariance_sqrt` uses, the order
/// of the symmetric eigensolver's backward error (Weyl perturbation bound), so
/// a computed eigenvalue within the floor of zero is indistinguishable from
/// zero. Its absolute term `N * denorm_min` keeps the floor from collapsing to
/// zero for a subnormal-scale M, below max|M_ij| = denorm_min / (floor_factor
/// eps); at and above that scale the floor is the relative term exactly. A
/// non-finite spectrum is not PSD.
///
/// @cite golub2013 -- Golub & Van Loan, "Matrix Computations", 4th ed., 2013, Sec. 8.1
template <typename Scalar, int N>
inline auto psd_decomposition(const Eigen::Matrix<Scalar, N, N>& M, Scalar floor_factor = Scalar{1}) -> psd_decomposition_result<Scalar, N>
{
    const Scalar floor = rounding_bound(N, M.cwiseAbs().maxCoeff(), floor_factor);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, N, N>> eig(symmetrize(M));
    const bool psd = (eig.eigenvalues().array() >= -floor).all();
    return {eig.eigenvalues(), eig.eigenvectors(), psd};
}

/// @brief Factor S with S*S^T equal to the decomposed matrix with its
/// rounding-level negative eigenvalues set to zero.
///
/// Unlike `covariance_sqrt`, a zero eigenvalue stays zero, which is what a
/// point-mass prior or a rank-one process noise requires.
template <typename Scalar, int N>
inline auto psd_factor(const psd_decomposition_result<Scalar, N>& d) -> Eigen::Matrix<Scalar, N, N>
{
    return d.eigenvectors * d.eigenvalues.cwiseMax(Scalar{0}).cwiseSqrt().asDiagonal();
}

/// @brief Solve s k = b for the transposed gain k of a one-output filter, with
/// s the 1 x 1 innovation covariance and b the 1 x N right-hand side.
///
/// This is bitwise Eigen 3.4's 1 x 1 column-pivoting QR solve, whose reflector
/// is the identity and whose triangular kernel divides one right-hand value
/// unless it is zero and multiplies several by 1 / s. Skipping the
/// factorization keeps gcc's optimized build from flagging its solve as an
/// over-read (-Wstringop-overread).
template <typename Derived>
inline auto one_output_gain_transpose(typename Derived::Scalar s, const Eigen::MatrixBase<Derived>& b) -> Eigen::Matrix<typename Derived::Scalar, 1, Derived::ColsAtCompileTime>
{
    using Scalar = typename Derived::Scalar;
    using row_t = Eigen::Matrix<Scalar, 1, Derived::ColsAtCompileTime>;
    const row_t rhs = b;
    if constexpr(Derived::ColsAtCompileTime == 1)
        return rhs(0) != Scalar{0} ? row_t{rhs / s} : rhs;
    else
        return rhs * (Scalar{1} / s);
}

}

#endif
