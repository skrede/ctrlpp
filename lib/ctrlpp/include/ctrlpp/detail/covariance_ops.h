#ifndef HPP_GUARD_CTRLPP_DETAIL_COVARIANCE_OPS_H
#define HPP_GUARD_CTRLPP_DETAIL_COVARIANCE_OPS_H

// @cite simon2006 -- Simon, "Optimal State Estimation", 2006, Sec. 5.5, p. 140
// @cite higham1988 -- Higham, "Computing a nearest symmetric positive semidefinite matrix", Linear Algebra Appl. 103, 1988
// @cite higham2002 -- Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed., 2002, Sec. 2.12, Eq. (2.8), and Sec. 3.1
// @cite golub2013 -- Golub & Van Loan, "Matrix Computations", 4th ed., 2013, Sec. 8.1 (symmetric eigenvalue perturbation)

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

#include <limits>
#include <cstdint>
#include <algorithm>

namespace ctrlpp::detail
{

// The bound `factor * count * eps * scale`, never below the underflow term `count * denorm_min`.
// `count` is in machine epsilons, eps = 2 u, derived by each caller from the convention that an
// n-term inner product carries at most gamma_n, about n u, of its scale (@cite higham2002, Sec. 3.1);
// `factor` is the order-one prefactor left to the caller. The absolute term is the gradual-underflow
// part of the standard model, fl(x op y) = (x op y)(1 + delta) + eta with |eta| <= denorm_min / 2
// (@cite higham2002, Sec. 2.12, Eq. (2.8)): without it the bound collapses to zero below
// scale = denorm_min / (factor eps), and at and above that scale the bound is the relative term exactly.
template <typename Scalar>
inline Scalar rounding_bound(std::int32_t count, Scalar scale, Scalar factor)
{
    const Scalar n = static_cast<Scalar>(count);
    return std::max(factor * n * std::numeric_limits<Scalar>::epsilon() * scale, n * std::numeric_limits<Scalar>::denorm_min());
}

// Symmetrizes M in place when its asymmetry is rounding drift, and refuses it, unchanged, when the
// asymmetry is structural. The count behind `factor = 1`: each triangle of a matrix formed by
// products carries at most N u of its scale, so the two triangles of one exactly symmetric matrix
// differ by at most 2 N u = N eps times it. Each differing pair is replaced in both positions by its
// midpoint M_ij + (M_ji - M_ij) / 2, which cannot overflow; an equal pair is never written, so a
// symmetric M is left bitwise untouched. A non-finite M is refused: an infinite entry would make the
// bound infinite too.
template <typename Scalar, int N>
inline bool symmetrize_within_drift(Eigen::Matrix<Scalar, N, N>& M, Scalar factor = Scalar{1})
{
    if(!M.allFinite())
        return false;
    const Scalar bound = rounding_bound(N, M.cwiseAbs().maxCoeff(), factor);
    if(!((M - M.transpose()).cwiseAbs().array() <= bound).all())
        return false;
    for(Eigen::Index j = 1; j < N; ++j)
        for(Eigen::Index i = 0; i < j; ++i)
            if(M(i, j) != M(j, i))
                M(i, j) = M(j, i) = M(i, j) + (M(j, i) - M(i, j)) / Scalar{2};
    return true;
}

// @cite simon2006 -- Simon, "Optimal State Estimation", 2006, Ch. 5 (numerical covariance conditioning)
template <typename Derived>
inline Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime> symmetrize(const Eigen::MatrixBase<Derived>& M)
{
    using Scalar = typename Derived::Scalar;
    return (Scalar{0.5} * (M + M.transpose())).eval();
}

template <typename Scalar, int N>
struct covariance_sqrt_result
{
    Eigen::Matrix<Scalar, N, N> factor;
    bool repaired;
};

// The unpivoted lower Cholesky factor squares back to P exactly; a pivoted LDLT's L*sqrt(D) would
// reconstruct a symmetric permutation of P instead. A P that is not positive definite is repaired
// to the nearest symmetric positive definite matrix (@cite higham1988) by clamping its eigenvalues
// up to `floor_factor * N * eps * max|P_ij|`, the order of the symmetric eigensolver's backward
// error (Weyl perturbation bound, @cite golub2013, Sec. 8.1), and `repaired` reports it.
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

template <typename Scalar, int N>
struct psd_decomposition_result
{
    Eigen::Matrix<Scalar, N, 1> eigenvalues;
    Eigen::Matrix<Scalar, N, N> eigenvectors;
    bool psd;
};

// Refuses rather than repairs. M is PSD when its smallest eigenvalue is at least
// -rounding_bound(N, max|M_ij|, floor_factor): the relative term is covariance_sqrt's floor, within
// which a computed eigenvalue is indistinguishable from zero, and the absolute term keeps the floor
// from collapsing to zero for a subnormal-scale M. A non-finite spectrum is not PSD.
// @cite golub2013 -- Golub & Van Loan, "Matrix Computations", 4th ed., 2013, Sec. 8.1
template <typename Scalar, int N>
inline psd_decomposition_result<Scalar, N> psd_decomposition(const Eigen::Matrix<Scalar, N, N>& M, Scalar floor_factor = Scalar{1})
{
    const Scalar floor = rounding_bound(N, M.cwiseAbs().maxCoeff(), floor_factor);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, N, N>> eig(symmetrize(M));
    const bool psd = (eig.eigenvalues().array() >= -floor).all();
    return {eig.eigenvalues(), eig.eigenvectors(), psd};
}

// Unlike covariance_sqrt, a zero eigenvalue stays zero, which is what a point-mass prior or a
// rank-one process noise requires.
template <typename Scalar, int N>
inline Eigen::Matrix<Scalar, N, N> psd_factor(const psd_decomposition_result<Scalar, N>& d)
{
    return d.eigenvectors * d.eigenvalues.cwiseMax(Scalar{0}).cwiseSqrt().asDiagonal();
}

// Bitwise Eigen 3.4's 1 x 1 column-pivoting QR solve of s k = b, whose reflector is the identity and
// whose triangular kernel divides one right-hand value unless it is zero and multiplies several by
// 1 / s. Skipping the factorization keeps gcc's optimized build from flagging its solve as an
// over-read (-Wstringop-overread).
template <typename Derived>
inline Eigen::Matrix<typename Derived::Scalar, 1, Derived::ColsAtCompileTime> one_output_gain_transpose(typename Derived::Scalar s, const Eigen::MatrixBase<Derived>& b)
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
