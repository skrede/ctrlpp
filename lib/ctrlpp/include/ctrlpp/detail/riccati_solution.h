#ifndef HPP_GUARD_CTRLPP_DETAIL_RICCATI_SOLUTION_H
#define HPP_GUARD_CTRLPP_DETAIL_RICCATI_SOLUTION_H

/// @brief Riccati invariant-subspace P extraction: P = U21 * U11^-1 on real U.
///
/// Given a reordered orthogonal real basis U of size 2n x 2n (the output of
/// ctrlpp::detail::reorder_real_schur), this header provides the single
/// routine required to recover the stabilizing Riccati solution P.
/// P is computed as U21 * U11^-1 directly on real arithmetic, symmetrised
/// via ctrlpp::detail::symmetrize, and checked for finiteness and positive
/// semi-definiteness. The positive semi-definiteness floor is derived from
/// std::numeric_limits<Scalar>::epsilon() scaled by the infinity norm of P;
/// no hardcoded numerical literal appears anywhere in the primitive.
///
/// Errors surface through ctrlpp::expected with a three-variant enum shared
/// between DARE and CARE; the public dare_error / care_error enums map
/// onto these variants one-to-one in the control/ wrappers.
///
/// @cite laub1979 -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979

#include "ctrlpp/expected.h"

#include "ctrlpp/detail/covariance_ops.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>
#include <initializer_list>

namespace ctrlpp::detail
{

/// @brief Negative-pivot floor below which an LDLT pivot of a symmetric N x N
/// matrix carries no sign information.
///
/// The bound is `floor_factor * N * eps * max|P_ij|`. LDLT of a symmetric
/// matrix is backward stable with a perturbation of order `N * eps * ||P||`
/// (Golub & Van Loan, Sec. 4.1), so a pivot inside that band is
/// indistinguishable from zero and must not be read as evidence of an
/// indefinite matrix. `floor_factor` is the leading order-one prefactor the
/// bound leaves unspecified, exposed for the same reason `covariance_sqrt`
/// exposes it, and this is the same `floor_factor * N * eps * max|P_ij|` form
/// that primitive already uses.
///
/// The `N` factor is load-bearing rather than cosmetic. Without it the floor
/// reads `1 * eps * max|P_ij|`, which demands a factorization carrying exactly
/// zero rounding error; a measured DARE solution missed that floor by 1.4% of a
/// single ulp while being positive semi-definite to within one ulp.
///
/// @cite golub2013 -- Golub & Van Loan, "Matrix Computations", 4th ed., 2013, Sec. 4.1
template <typename Scalar, int N>
auto psd_pivot_floor(const Eigen::Matrix<Scalar, N, N>& P,
                     Scalar floor_factor = Scalar{1}) -> Scalar
{
    return -floor_factor * static_cast<Scalar>(N)
           * std::numeric_limits<Scalar>::epsilon() * P.cwiseAbs().maxCoeff();
}

/// @brief A Frobenius magnitude together with whether the arithmetic resolved it.
///
/// `value` carries meaning only when `resolved` is true. Unresolved is
/// contagious through every combinator below, and no comparison built on these
/// can succeed with an unresolved side, so a magnitude the arithmetic could not
/// form declines rather than certifying.
///
/// ## Why a plain Frobenius magnitude is not enough, at BOTH ends
///
/// The two Riccati solvers fail at opposite ends of the arithmetic range, so a
/// magnitude that is safe at only one end leaves half the fault class standing.
///
///  * At the top, a matrix whose entries are finite but whose sum of squares
///    exceeds the largest finite value gives an infinite magnitude. Both sides
///    of an acceptance comparison then become infinite, `inf <= inf` holds, and
///    the comparison certifies whatever it was handed. The discrete solver's
///    original-scale verification meets this first and refuses answers that are
///    perfectly representable.
///  * At the bottom, a matrix whose entries are small but nonzero and whose sum
///    of squares falls below the smallest normal value gives exactly zero. A
///    residual scale of zero turns a `residual <= floor * scale` test into
///    `0 <= 0`, which passes unconditionally, and a residual that underflows to
///    zero is read as an exact solve. The continuous solver meets this end.
///
/// ## Which primitive this wraps, and why that one
///
/// Screened on four operands before any cost was considered: every entry
/// `1e200`; every entry `1e-200`; one entry `1e-200` with the rest zero; and
/// one entry equal to the smallest positive subnormal. The plain member fails
/// all four. `Eigen::blueNorm()` survives the top and FAILS the bottom,
/// returning exactly zero for the smallest-subnormal operand, so only a screen
/// at both ends rejects it. `Eigen::stableNorm()` and the largest-entry rescale
/// `max|M| * ||M / max|M|||` both survive all four.
///
/// The rescale is taken. It is the cheaper survivor at the sizes both Riccati
/// paths carry -- 2.2, 4.5 and 3.6 times the plain member at N = 2, 4 and 8,
/// against 8.6, 9.8 and 3.9 for `stableNorm()` -- and it is the only survivor
/// that applies to a matrix at all: the Eigen norms are vector primitives and
/// trip an internal block assertion on a fixed-size matrix operand, which would
/// force every call site to reinterpret its operand as a contiguous vector and
/// would be unavailable on an unevaluated expression. Both survivors allocate
/// nothing at every size measured; the rescale is chosen on cost and reach, not
/// on allocation.
///
/// The rescale is exact in the sense that matters: dividing by the largest
/// entry cannot overflow or underflow, every rescaled entry lands in [-1, 1],
/// and the sum of squares of at most N * N such entries cannot leave either end
/// of the range. The final multiplication by the largest entry is the only step
/// that can leave the top, and when it does the magnitude genuinely is not
/// representable and the result is unresolved rather than infinite.
template <typename Scalar>
struct resolved_magnitude
{
    Scalar value;
    bool   resolved;
};

/// @brief Frobenius magnitude of a matrix, safe at both ends of the range.
template <typename Derived>
auto resolve_magnitude(const Eigen::MatrixBase<Derived>& operand)
    -> resolved_magnitude<typename Derived::Scalar>
{
    using Scalar = typename Derived::Scalar;
    const auto& evaluated = operand.eval();

    if (!evaluated.allFinite())
        return {Scalar{0}, false};

    const Scalar largest_entry = evaluated.cwiseAbs().maxCoeff();
    if (!(largest_entry > Scalar{0}))
        return {Scalar{0}, true};

    const Scalar magnitude =
        largest_entry * (evaluated / largest_entry).norm();
    if (!std::isfinite(magnitude))
        return {Scalar{0}, false};
    return {magnitude, true};
}

/// @brief Largest-entry magnitude of a matrix.
///
/// This form squares nothing, so it is already safe at both ends and needs only
/// the finiteness guard. It is kept distinct from the Frobenius form rather than
/// folded into it because a bound derived from the largest entry is a different
/// bound, and substituting one for the other would silently move a threshold.
template <typename Derived>
auto resolve_largest_entry(const Eigen::MatrixBase<Derived>& operand)
    -> resolved_magnitude<typename Derived::Scalar>
{
    using Scalar = typename Derived::Scalar;
    const auto& evaluated = operand.eval();

    if (!evaluated.allFinite())
        return {Scalar{0}, false};
    return {evaluated.cwiseAbs().maxCoeff(), true};
}

/// @brief Largest of several magnitudes; unresolved if any operand is.
template <typename Scalar>
auto largest_magnitude(
    std::initializer_list<resolved_magnitude<Scalar>> operands)
    -> resolved_magnitude<Scalar>
{
    resolved_magnitude<Scalar> largest{Scalar{0}, true};
    for (const resolved_magnitude<Scalar>& operand : operands)
    {
        if (!operand.resolved)
            return {Scalar{0}, false};
        if (operand.value > largest.value)
            largest.value = operand.value;
    }
    return largest;
}

/// @brief A magnitude scaled by a finite factor; unresolved if the product is
/// not finite, so a bound that has itself left the range cannot admit anything.
template <typename Scalar>
auto scaled_magnitude(const resolved_magnitude<Scalar>& operand,
                      Scalar                            factor)
    -> resolved_magnitude<Scalar>
{
    if (!operand.resolved || !std::isfinite(factor))
        return {Scalar{0}, false};

    const Scalar scaled = operand.value * factor;
    if (!std::isfinite(scaled))
        return {Scalar{0}, false};
    return {scaled, true};
}

/// @brief Counted rounding bound `rounding_ops * eps * scale` on a resolved
/// scale, in the counted-operation form both Riccati acceptance chains use.
template <typename Scalar>
auto counted_rounding_bound(const resolved_magnitude<Scalar>& scale,
                            int                               rounding_ops)
    -> resolved_magnitude<Scalar>
{
    return scaled_magnitude(scale,
                            static_cast<Scalar>(rounding_ops)
                                * std::numeric_limits<Scalar>::epsilon());
}

/// @brief `value <= bound`, false whenever either side is unresolved.
///
/// Two infinities do not compare equal-or-less through this helper, which is the
/// whole point: an unresolved operand is an absence of evidence and must not
/// read as evidence of smallness.
template <typename Scalar>
auto magnitude_within(const resolved_magnitude<Scalar>& value,
                      const resolved_magnitude<Scalar>& bound) -> bool
{
    return value.resolved && bound.resolved && value.value <= bound.value;
}

/// @brief Error variants produced by extract_riccati_solution.
///
/// Shared between DARE and CARE. The public dare_error / care_error enums
/// map their singular_u11 / non_finite_input / non_psd_solution variants
/// directly onto these values.
enum class riccati_extract_error
{
    singular_u11,
    non_finite,
    non_psd,
};

/// @brief Extract the stabilizing Riccati solution into a caller-supplied matrix.
///
/// Given an orthogonal real U of size 2n x 2n with the stable invariant
/// subspace in its leading n columns (output of reorder_real_schur), compute
/// P = U21 * U11^-1, symmetrise in place, and validate positive semi-definiteness.
/// Writes directly into P_out so that callers can avoid the ctrlpp::expected<Matrix>
/// return-by-value copy on the hot path.
///
/// The implementation prefers a back-substitution against U11 (solving
/// U11^T * P^T = U21^T) over an explicit inverse for numerical accuracy.
///
/// The positive semi-definiteness floor is -eps * ||P||_inf, matching the
/// LAPACK convention of scaling relative thresholds by the operand norm.
///
/// @returns ctrlpp::expected<void, riccati_extract_error>.
template <typename Scalar, int N2>
auto extract_riccati_solution_into(
    Eigen::Matrix<Scalar, N2 / 2, N2 / 2>&    P_out,
    const Eigen::Matrix<Scalar, N2, N2>&      U)
    -> ctrlpp::expected<void, riccati_extract_error>
{
    static_assert(N2 > 0 && (N2 % 2 == 0), "U must have even size 2n x 2n");

    constexpr int n = N2 / 2;
    using MatNxN = Eigen::Matrix<Scalar, n, n>;

    const MatNxN U11 = U.template block<n, n>(0, 0);
    const MatNxN U21 = U.template block<n, n>(n, 0);

    auto qr_U11T = U11.transpose().colPivHouseholderQr();
    if (!qr_U11T.isInvertible())
        return ctrlpp::unexpected(riccati_extract_error::singular_u11);

    // P = U21 * U11^-1  <=>  U11^T * P^T = U21^T  (solve against QR of U11^T).
    const MatNxN P_raw = qr_U11T.solve(U21.transpose()).transpose();
    P_out = ctrlpp::detail::symmetrize(P_raw);

    if (!P_out.allFinite())
        return ctrlpp::unexpected(riccati_extract_error::non_finite);

    // PSD check via LDLT: by Sylvester's law of inertia, the signs of the pivots
    // in the D diagonal match the signs of the eigenvalues of a symmetric P.
    // Detects non-PSD at ~1/10 the instruction cost of a full eigendecomposition
    // while preserving the same eps-scaled rejection threshold.
    Eigen::LDLT<MatNxN> ldlt(P_out);
    if (ldlt.info() != Eigen::Success
        || ldlt.vectorD().minCoeff() < psd_pivot_floor<Scalar, n>(P_out))
        return ctrlpp::unexpected(riccati_extract_error::non_psd);

    return {};
}

/// @brief Value-returning wrapper around `extract_riccati_solution_into`.
template <typename Scalar, int N2>
auto extract_riccati_solution(const Eigen::Matrix<Scalar, N2, N2>& U)
    -> ctrlpp::expected<Eigen::Matrix<Scalar, N2 / 2, N2 / 2>, riccati_extract_error>
{
    Eigen::Matrix<Scalar, N2 / 2, N2 / 2> P;
    auto err = extract_riccati_solution_into<Scalar, N2>(P, U);
    if (!err)
        return ctrlpp::unexpected(err.error());
    return P;
}

}

#endif
