#ifndef HPP_GUARD_CTRLPP_DETAIL_COVARIANCE_OPS_H
#define HPP_GUARD_CTRLPP_DETAIL_COVARIANCE_OPS_H

/// @brief Covariance-matrix utilities (symmetrisation, matrix square root).
///
/// @cite simon2006 -- Simon, "Optimal State Estimation", 2006, Ch. 5 (covariance conditioning / Joseph form)
/// @cite higham1988 -- Higham, "Computing a nearest symmetric positive semidefinite matrix", Linear Algebra Appl. 103, 1988
/// @cite golub2013 -- Golub & Van Loan, "Matrix Computations", 4th ed., 2013, Sec. 8.1 (symmetric eigenvalue perturbation)

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

#include <limits>

namespace ctrlpp::detail
{

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

}

#endif
