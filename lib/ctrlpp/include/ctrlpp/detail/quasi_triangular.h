#ifndef HPP_GUARD_CTRLPP_DETAIL_QUASI_TRIANGULAR_H
#define HPP_GUARD_CTRLPP_DETAIL_QUASI_TRIANGULAR_H

/// @brief Block structure and eigenvalues of a real quasi-triangular factor,
/// defined once.
///
/// A real Schur factor is block upper triangular with 1x1 and 2x2 diagonal
/// blocks. Two questions are asked of it repeatedly -- where does a block start
/// and end, and what are its eigenvalues -- and both have a trap in them that a
/// second copy is likely to fall into.
///
/// The first trap is the diagonal. A 2x2 block's eigenvalues are the half-trace
/// plus and minus the square root of the discriminant, NOT the two diagonal
/// entries: for a complex-conjugate pair the diagonal entries are not the shared
/// real part unless the block has been standardized, and for a real pair they are
/// neither root. Reading the diagonal is wrong in both directions, and which
/// direction it is wrong in decides whether a caller over-refuses or certifies
/// something it should not.
///
/// The second trap is the subdiagonal entry. It is what separates a 2x2 block
/// from two 1x1 blocks, and on a COMPUTED factor it is a rounded quantity, so
/// testing it against exact zero asks a question the arithmetic cannot answer.
/// The test here is the LAPACK one: significant when it exceeds unit roundoff
/// times the larger of the whole factor's magnitude and the block's own diagonal
/// magnitude. The global term keeps a block inside a badly scaled factor from
/// being split on rounding noise; the local term keeps a small block inside a
/// large factor from being fused when its own entries say otherwise.
///
/// Callers that already know a block's eigenvalues are real or complex do not
/// re-derive it: `block_spectrum::real_pair` carries the discriminant's verdict.
///
/// @cite golub2013 -- Golub & Van Loan, "Matrix Computations", 4th ed., 2013, Sec. 7.4.1
/// @cite bai_demmel_1993 -- Bai & Demmel, "On swapping diagonal blocks in real Schur form", 1993

#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <complex>
#include <algorithm>

namespace ctrlpp::detail
{

/// @brief The two eigenvalues of one diagonal block of a real Schur factor.
///
/// For a 1x1 block both entries hold the same real eigenvalue. For a 2x2 block
/// `first` is the root with the larger real part, so a caller testing a
/// half-plane membership tests `first` and needs no branch on `real_pair`.
template <typename Scalar>
struct block_spectrum
{
    std::complex<Scalar> first;
    std::complex<Scalar> second;
    bool                 real_pair;
};

/// @brief Size of the diagonal block starting at `pos`: 1 or 2.
///
/// `factor_scale` is the largest absolute entry of the whole factor. Passing it
/// in rather than recomputing it keeps a walk over all blocks linear in the
/// factor's size rather than quadratic.
template <typename Scalar, int N>
auto quasi_triangular_block_size(const Eigen::Matrix<Scalar, N, N>& T,
                                 int                                pos,
                                 Scalar factor_scale) -> int
{
    using std::abs;

    if (pos + 1 >= N)
        return 1;

    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    const Scalar local_scale =
        std::max(abs(T(pos, pos)), abs(T(pos + 1, pos + 1)));
    return (abs(T(pos + 1, pos)) > eps * std::max(factor_scale, local_scale))
               ? 2
               : 1;
}

/// @brief Eigenvalues of the diagonal block at `pos`, from trace and
/// determinant.
template <typename Scalar, int N>
auto quasi_triangular_block_spectrum(const Eigen::Matrix<Scalar, N, N>& T,
                                     int                                pos,
                                     int block_size) -> block_spectrum<Scalar>
{
    using std::sqrt;

    if (block_size == 1)
    {
        const std::complex<Scalar> eigenvalue(T(pos, pos), Scalar{0});
        return {eigenvalue, eigenvalue, true};
    }

    const Scalar a = T(pos,     pos);
    const Scalar b = T(pos,     pos + 1);
    const Scalar c = T(pos + 1, pos);
    const Scalar d = T(pos + 1, pos + 1);
    const Scalar half_trace   = (a + d) / Scalar{2};
    const Scalar determinant  = a * d - b * c;
    const Scalar discriminant = half_trace * half_trace - determinant;

    if (discriminant >= Scalar{0})
    {
        const Scalar root = sqrt(discriminant);
        return {std::complex<Scalar>(half_trace + root, Scalar{0}),
                std::complex<Scalar>(half_trace - root, Scalar{0}),
                true};
    }

    const Scalar imaginary_part = sqrt(-discriminant);
    return {std::complex<Scalar>(half_trace,  imaginary_part),
            std::complex<Scalar>(half_trace, -imaginary_part),
            false};
}

/// @brief Backward-error margin on an eigenvalue's real part read from a real
/// Schur factor: the factored matrix's dimension times unit roundoff times the
/// factor's largest entry.
///
/// A backward-stable factorization returns the exact factor of a nearby matrix,
/// and an eigenvalue of an N-by-N factor carries a perturbation on the order of
/// N times unit roundoff times the factor's magnitude. An open-half-plane
/// predicate that does not subtract this cannot distinguish an eigenvalue on the
/// axis from one the arithmetic placed there.
///
/// Every left-half-plane predicate in the continuous solver uses this one
/// derivation: the plain-Schur and balanced-Schur reordering predicates on the
/// 2n-by-2n Hamiltonian factor, and the acceptance rule on the n-by-n
/// closed-loop factor. The acceptance rule reaches it through the
/// resolved-magnitude form, which additionally declines when the factor's
/// magnitude cannot be resolved, but the counted expression is the same.
template <typename Scalar, int N>
auto schur_eigenvalue_margin(const Eigen::Matrix<Scalar, N, N>& T) -> Scalar
{
    return Scalar{N} * std::numeric_limits<Scalar>::epsilon()
           * T.cwiseAbs().maxCoeff();
}

/// @brief True when every eigenvalue of the factor has real part below
/// `-margin`.
///
/// The binding eigenvalue of a block is the one with the larger real part, and
/// `block_spectrum::first` is that one in both the real-pair and the
/// conjugate-pair case, so the test carries no branch on which case it is in.
/// That is the whole point of asking the discriminant: a 2x2 block whose
/// eigenvalues are two REAL roots of opposite sign has a half-trace strictly
/// left of the margin while its larger root sits in the open right half-plane,
/// and a walk that reads the half-trace certifies it.
template <typename Scalar, int N>
auto quasi_triangular_spectrum_strictly_left_of(
    const Eigen::Matrix<Scalar, N, N>& T,
    Scalar                             margin) -> bool
{
    const Scalar factor_scale = T.cwiseAbs().maxCoeff();

    int index = 0;
    while (index < N)
    {
        const int block_size =
            quasi_triangular_block_size<Scalar, N>(T, index, factor_scale);
        const block_spectrum<Scalar> block =
            quasi_triangular_block_spectrum<Scalar, N>(T, index, block_size);
        if (!(block.first.real() < -margin))
            return false;
        index += block_size;
    }
    return true;
}

}

#endif
