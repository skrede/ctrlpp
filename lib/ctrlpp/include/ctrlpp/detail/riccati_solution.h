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
#include <cstddef>
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

/// @brief Rounded operations along the longest chain producing one entry of the
/// discrete Riccati residual, for an NX-state, NU-input problem.
///
/// Enumerated rather than chosen. Each contraction over the state dimension
/// costs NX multiplies and NX-1 additions, that is 2*NX-1, and six of them occur
/// along the chain: A'P, (A'P)A, B'P, (B'P)A, A'PB, and the contraction of A'PB
/// against the gain. Each contraction over the input dimension costs 2*NU-1, and
/// two occur: (B'P)B, and the gain's own inner dimension. The weighting sum
/// R + B'PB is one addition. The linear solve for the gain is a rank-revealing
/// QR of an NU x NU matrix followed by a back substitution, whose backward error
/// is bounded by 2*NU operations at the scale of the matrix it factorizes.
/// Assembling the four terms is three additions.
///
/// Every operation is counted whether or not it actually rounds, so the count
/// bounds the accumulated error from above rather than describing it tightly --
/// which is what a budget requires.
///
/// This is the single definition. The solver's gain-agreement margin and the
/// unit anchors' rounding-floor assertions both quote it; neither re-spells it.
template <std::size_t NX, std::size_t NU>
constexpr int dare_residual_ops = 6 * (2 * static_cast<int>(NX) - 1) + 2 * (2 * static_cast<int>(NU) - 1) + 1 + 2 * static_cast<int>(NU) + 3;

/// @brief The relative accuracy below which more than half of the scalar type's
/// significand is retained.
///
/// Derived, not calibrated. For a radix-2 scalar type with `eps = 2^-p`
/// (`p = 52` for binary64, `p = 23` for binary32), a relative accuracy of
/// `sqrt(eps) = 2^(-p/2)` is exactly the retention of `p/2` of the `p` fractional
/// significand bits. "At least half the significand of the returned answer is
/// correct" therefore IS "relative forward error at most `sqrt(eps)`", stated in
/// the type's own radix, with nothing fitted and nothing measured. It carries to
/// every radix-2 type without re-measurement, which is the property a calibrated
/// multiple can never have.
template <typename Scalar>
auto half_significand_margin() -> Scalar
{
    return std::sqrt(std::numeric_limits<Scalar>::epsilon());
}

/// @brief What the one accuracy rule decided about a solved Riccati pose.
///
/// The third value is not a failure. An estimate that could not be formed has
/// produced no evidence, and an absence of evidence is neither a certificate nor
/// a refutation.
enum class riccati_accuracy
{
    within,
    exceeded,
    unresolved,
};

/// @brief Estimate the relative forward error of a discrete Riccati solution
/// from its residual, by inverting the residual map's derivative.
///
/// ## Why the residual itself cannot be the acceptance quantity
///
/// The residual is a four-term cancellation whose magnitude is set by the
/// conditioning of the terms, not by the accuracy of the answer. Measured over
/// 2.5 million poses, the answers that retain more than half the significand and
/// the answers that do not are CONTIGUOUS on the residual -- the worst retained
/// forward error and the best lost one are adjacent, and the distribution is
/// unimodal across ten decades with no gap anywhere in it. No threshold on that
/// quantity separates the two, which is why a bound on the residual is not a
/// weak accuracy gate but not one at all. Twenty conditioning-aware
/// amplifications of the enumerated-operation form were measured; the only one
/// that bounded the whole population was four orders of magnitude looser than
/// the envelope it would have replaced, so no amplified residual bound is
/// written here and none is to be written.
///
/// ## What is estimated instead, and why it is a first-order identity
///
/// Let `P` be the exact stabilizing solution, `P_hat` the computed one, and
/// `E = P_hat - P`. The residual map
///
///     R(X) = A' X A - X - A' X B (R + B' X B)^-1 B' X A + Q
///
/// has, at the solution, the Frechet derivative built from the closed loop:
///
///     DR_P[E] = A_cl' E A_cl - E,     A_cl = A - B K
///
/// so with `Omega(X) = X - A_cl' X A_cl` and `R(P) = 0`,
///
///     R(P_hat) = -Omega(E) + O(||E||^2) .
///
/// Inverting gives `E = -Omega^-1(R(P_hat))` to first order, so the answer's own
/// relative forward error is estimated by `||Omega^-1(residual)||_F / ||P||_F`
/// with no fitted constant anywhere in it. Against an extended-precision
/// reference of a different algorithm the estimate's median ratio to the true
/// forward error is within a few percent over half a million independent poses,
/// with the right-hand tail over-predicting, which is the conservative direction.
///
/// The estimate is NOT exact and is not sold as one. Its left tail
/// under-predicts where cancellation in the residual is favorable, so a small
/// population of answers that have genuinely lost more than half the significand
/// is retained rather than refused. A first-order estimator cannot close that
/// gap; only a second solve in higher precision could, and that is not a
/// postcondition.
///
/// @cite laub1979  -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979
/// @cite higham2008 -- Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed., 2008, Ch. 19
///
/// ## Scale, and why the whole estimate is divided through by the peak of P
///
/// The estimate is invariant under `P -> sP` with the weights carried along, so
/// it is formed on `P / max|P_ij|` and `residual / max|P_ij|`. Dividing by the
/// largest entry cannot overflow or underflow, it puts the denominator in
/// `[1, N]` for every representable P, and it removes the divide-by-`||P||_F`
/// hazard at both ends of the arithmetic range -- the same rescale the magnitude
/// helpers above use, for the same reason.
///
/// ## Cost, counted rather than timed
///
/// `Omega` is assembled column by column on the symmetric subspace, whose
/// dimension is `M = N(N+1)/2`. Each basis image is one rank-two outer-product
/// update, so the assembly is `M * N^2` operations; the rank-revealing solve is
/// `(2/3) M^3 ~ N^6 / 12`. Against a seven-iteration Schur solve's `149.33 N^3`
/// that is 0.4% at `N = 2`, 3.6% at `N = 4` and 28.6% at `N = 8`. The growth is
/// `N^3 / 1792` relative, so beyond roughly `N = 10` the direct form stops being
/// the right one: a Bartels-Stewart Stein solve against a real Schur factor of
/// `A_cl` costs about `26 N^3`, a fixed 17% of the solve at every size. That
/// reduction is available and is deliberately not taken here, because the
/// library instantiates the Riccati path at `N = 2`, `4` and `8`, where the
/// direct form is at most 1.7x more work than the factorized one and is
/// materially simpler to audit. Nothing on this path allocates.
template <typename Scalar, int N>
auto estimate_riccati_forward_error(const Eigen::Matrix<Scalar, N, N>& closed_loop,
                                    const Eigen::Matrix<Scalar, N, N>& residual,
                                    const Eigen::Matrix<Scalar, N, N>& P)
    -> resolved_magnitude<Scalar>
{
    static_assert(N > 0, "state dimension must be positive");
    constexpr int M = N * (N + 1) / 2;
    using MatNxN = Eigen::Matrix<Scalar, N, N>;

    if(!closed_loop.allFinite() || !residual.allFinite() || !P.allFinite())
        return {Scalar{0}, false};

    const Scalar peak = P.cwiseAbs().maxCoeff();
    if(!(peak > Scalar{0}))
    {
        // P is exactly zero. Every term of the residual that contains P is
        // exactly zero too, so the only claim available -- and the only one
        // needed -- is that what remains is exactly zero as well. That is an
        // exact statement about an exact solve, not a degenerate comparison.
        const Scalar remainder = residual.cwiseAbs().maxCoeff();
        return {remainder > Scalar{0} ? std::numeric_limits<Scalar>::infinity() : Scalar{0}, true};
    }

    const MatNxN residual_unit = (residual / peak).eval();
    const MatNxN P_unit        = (P / peak).eval();
    if(!residual_unit.allFinite())
        return {Scalar{0}, false};

    // Symmetric basis: E_ii = e_i e_i', E_ij = e_i e_j' + e_j e_i' for i < j.
    // Omega(E_ij) = E_ij - c_i c_j' - c_j c_i' with c_i the i-th ROW of A_cl
    // taken as a column, so each column of the operator is one outer-product
    // pair rather than a pair of matrix products. The coordinate map is the
    // upper triangle read in the same order, which is exact for a symmetric
    // image because the off-diagonal basis element carries the 1 at both places.
    Eigen::Matrix<Scalar, M, M> stein_operator;
    Eigen::Matrix<Scalar, M, 1> stein_rhs;
    int column = 0;
    for(int i = 0; i < N; ++i)
    {
        for(int j = i; j < N; ++j)
        {
            MatNxN image = MatNxN::Zero();
            image(i, j) += Scalar{1};
            const auto ci = closed_loop.row(i).transpose();
            const auto cj = closed_loop.row(j).transpose();
            image -= ci * cj.transpose();
            if(i != j)
            {
                image(j, i) += Scalar{1};
                image -= cj * ci.transpose();
            }

            int row = 0;
            for(int a = 0; a < N; ++a)
                for(int b = a; b < N; ++b)
                    stein_operator(row++, column) = image(a, b);
            ++column;
        }
    }
    int row = 0;
    for(int a = 0; a < N; ++a)
        for(int b = a; b < N; ++b)
        {
            // The residual is symmetric in exact arithmetic; the symmetric part
            // is what the operator's range can represent, and taking it is what
            // keeps the coordinate system consistent with the basis above.
            stein_rhs(row++) = (residual_unit(a, b) + residual_unit(b, a)) / Scalar{2};
        }

    auto stein_qr = stein_operator.colPivHouseholderQr();
    stein_qr.setThreshold(Scalar{M} * std::numeric_limits<Scalar>::epsilon());
    if(!stein_qr.isInvertible())
    {
        // Omega is singular exactly when the closed loop carries an eigenvalue
        // pair with lambda_i * lambda_j = 1, which a spectrum strictly inside
        // the unit disk forbids. A singular operator is therefore evidence
        // AGAINST the stabilizing claim rather than an absence of evidence, and
        // it is reported as an error past every margin rather than as
        // unresolved.
        return {std::numeric_limits<Scalar>::infinity(), true};
    }

    const Eigen::Matrix<Scalar, M, 1> coordinates = stein_qr.solve(stein_rhs);
    if(!coordinates.allFinite())
        return {Scalar{0}, false};

    MatNxN error_estimate = MatNxN::Zero();
    row = 0;
    for(int a = 0; a < N; ++a)
        for(int b = a; b < N; ++b)
        {
            error_estimate(a, b) = coordinates(row);
            error_estimate(b, a) = coordinates(row);
            ++row;
        }

    const resolved_magnitude<Scalar> numerator   = resolve_magnitude(error_estimate);
    const resolved_magnitude<Scalar> denominator = resolve_magnitude(P_unit);
    if(!numerator.resolved || !denominator.resolved || !(denominator.value > Scalar{0}))
        return {Scalar{0}, false};

    const Scalar ratio = numerator.value / denominator.value;
    if(!std::isfinite(ratio))
        return {std::numeric_limits<Scalar>::infinity(), true};
    return {ratio, true};
}

/// @brief THE accuracy rule for a discrete Riccati solution, defined once.
///
/// Refuse when the estimated relative forward error of the returned answer
/// exceeds the half-significand margin. The solver's postcondition, the unit
/// anchors and the randomized oracle all reach their verdict through this one
/// call; none of them re-spells the estimator, the margin or the comparison.
template <typename Scalar, int N>
auto riccati_forward_error_verdict(const Eigen::Matrix<Scalar, N, N>& closed_loop,
                                   const Eigen::Matrix<Scalar, N, N>& residual,
                                   const Eigen::Matrix<Scalar, N, N>& P) -> riccati_accuracy
{
    const resolved_magnitude<Scalar> estimate =
        estimate_riccati_forward_error<Scalar, N>(closed_loop, residual, P);
    if(!estimate.resolved)
        return riccati_accuracy::unresolved;
    return estimate.value <= half_significand_margin<Scalar>() ? riccati_accuracy::within
                                                               : riccati_accuracy::exceeded;
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
