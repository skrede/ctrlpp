#ifndef HPP_GUARD_CTRLPP_DETAIL_RICCATI_SOLUTION_H
#define HPP_GUARD_CTRLPP_DETAIL_RICCATI_SOLUTION_H

/// @brief Riccati invariant-subspace P extraction: P = U21 * U11^-1 on real U.
///
/// Given a reordered orthogonal real basis U of size 2n x 2n (the output of
/// ctrlpp::detail::reorder_real_schur), this header provides the single
/// routine required to recover the stabilizing Riccati solution P.
/// P is computed as U21 * U11^-1 directly on real arithmetic, symmetrized
/// via ctrlpp::detail::symmetrize, and checked for finiteness and positive
/// semi-definiteness. The positive semi-definiteness floor is the bound
/// ctrlpp::detail::psd_pivot_floor computes: the floor factor times the state
/// dimension times std::numeric_limits<Scalar>::epsilon() times the LARGEST
/// ABSOLUTE ENTRY of P, negated. It is neither a row-sum norm nor
/// dimension-free; that function's docblock carries the derivation and the
/// measured solution the dimension factor exists for, and is not restated here.
/// No calibrated constant appears anywhere in the primitive: every numeric in
/// it is either a structural count fixed by the problem's shape or the
/// order-one prefactor the backward-error bound leaves unspecified, which is
/// exposed as a defaulted parameter rather than fixed in the body.
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
/// ## Arithmetic, counted rather than timed
///
/// `Omega` is assembled column by column on the symmetric subspace, whose
/// dimension is `M = N(N+1)/2`. Each basis image is one rank-two outer-product
/// update. Enumerating them gives `N` diagonal columns at `N^2` each plus
/// `N(N-1)/2` off-diagonal columns at `2 N^2` each, which is exactly `N^4` --
/// an off-diagonal column carries a SECOND outer product, so the shorter
/// `M * N^2` reading of the same loop undercounts it by a factor of two. The
/// rank-revealing solve is `(2/3) M^3`, and the coordinate round-trip and
/// back-substitution are `O(M^2)`. Against a seven-iteration Schur solve's
/// `149.33 N^3`, the whole estimate is 3.6% of the solve at `N = 2`, 10.7% at
/// `N = 4`, 24.5% at `N = 6` and 47.7% at `N = 8`. Nothing on this path
/// allocates.
///
/// A Bartels-Stewart Stein solve against a real Schur factor of `A_cl` needs
/// only `O(N^2)` storage and, counted the same way term by term, about
/// `18 N^3` operations: cheaper than the form here by 2.0x at `N = 6` and 3.9x
/// at `N = 8`, and 3.3x MORE expensive at `N = 2`, where the crossover has not
/// yet happened. It is not taken, and the reason is accuracy rather than cost.
/// Measured against an extended-precision reference of a different algorithm on
/// the same population this gate is justified by, the factorized form refused
/// four fewer wrong answers and its worst escaped forward error was 4.9x
/// larger, because the solvability check it can carry -- a pivot ratio local to
/// each `p*q <= 4` sub-system -- fires on a set disjoint from the one the
/// rank-revealing decomposition's own criterion fires on. Its accuracy is also
/// unmeasured at exactly the dimensions where its arithmetic advantage appears.
/// The reduction is real and the evidence for it is not, so what ships is the
/// form whose accuracy is already established.
///
/// ## Storage and stack, measured rather than asserted
///
/// STORAGE GROWS AS THE FOURTH POWER OF THE STATE DIMENSION, and the in-place
/// decomposition below does not change that. The operator is `M x M`, which is
/// asymptotically `N^4 / 4` entries; factorizing into its own storage removes
/// the second live copy of it, not the exponent. This is a bounded improvement
/// with a stated horizon, not an asymptotic answer.
///
/// Measured with `-fstack-usage` at `-O2` on x86-64, this function's own frame
/// is 864, 2,400, 6,432 and 15,056 bytes at `N = 2, 4, 6, 8`, against 896,
/// 3,168, 9,936 and 25,392 for the copying form it replaces -- a 1.69x
/// reduction at `N = 8`.
///
/// What a hard-real-time caller budgets is not this frame but the peak of the
/// WHOLE discrete solve chain, which `-fstack-usage` cannot give because it
/// attributes nothing to callees. Measured directly, by painting a region below
/// the frame and reading back the deepest disturbed word (harness floor zero on
/// every configuration), that peak is 5,352, 11,688, 23,144 and 42,760 bytes at
/// `N = 2, 4, 6, 8` for an `N`-state, 3-input pose. So, strictly and with no
/// margin left for the caller's own frames, interrupt context or RTOS overhead:
///
///     task stack   supported maximum N
///      4 KB        none, not even N = 2
///      8 KB        N <= 2
///     16 KB        N <= 4
///     32 KB        N <= 6
///     48 KB        N <= 8
///     64 KB        N <= 8
///
/// THE ESTIMATOR IS NOT WHAT DECIDES THE SMALL-STACK ANSWER. The same chain with
/// no accuracy estimate on it at all still peaks at 17,112 bytes at `N = 6` and
/// 28,168 at `N = 8`, so on a 4-16 KB task stack the supported maximum is
/// `N = 4` whether this estimate is formed or not. Removing it entirely would
/// buy no additional configuration below 32 KB. The construction here changes
/// the supported maximum in exactly one band, at 48 KB.
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

    // The decomposition is declared over a reference to the operator's own
    // storage and therefore factorizes IN PLACE, destroying `stein_operator`.
    // This is the only form in which Eigen writes the factorization into the
    // caller's array; both `stein_operator.colPivHouseholderQr()` and a
    // `ColPivHouseholderQR<Matrix<M, M>>` with a separate `compute()` hold an
    // `M x M` member of their own and leave the operand live beside it, so
    // either of those keeps two `M x M` arrays alive at once. The operator is
    // not read after this point, so there is nothing to preserve. The
    // arithmetic is unchanged: the same Householder sequence runs over the same
    // entries, and the solve is bit-for-bit what the copying form produces.
    Eigen::ColPivHouseholderQR<Eigen::Ref<Eigen::Matrix<Scalar, M, M>>> stein_qr(stein_operator);
    stein_qr.setThreshold(Scalar{M} * std::numeric_limits<Scalar>::epsilon());
    if(!stein_qr.isInvertible())
    {
        // AN OPERATOR TOO ILL-CONDITIONED TO INVERT IS AN ABSENCE OF EVIDENCE.
        //
        // Omega is EXACTLY singular only when the closed loop carries an
        // eigenvalue pair with lambda_i * lambda_j = 1, which a spectrum
        // strictly inside the unit disk forbids. The test above is not that
        // test. `isInvertible()` at a RELATIVE threshold reports rank deficiency
        // whenever the operator's condition number exceeds the reciprocal of
        // `M * eps`, which for Omega(X) = X - A_cl' X A_cl means roughly
        // 1 - |lambda|^2 < M * eps. The stabilizing claim is established
        // elsewhere against a margin of N * eps * ||A_cl||, so there is a band
        // -- for binary64 at N = 8, |lambda| between 1 - 4e-15 and 1 - 4.4e-16
        // -- in which the spectrum check passes, the answer may be perfectly
        // good, and this decomposition still cannot form the estimate. A
        // conditioning statement about the arithmetic is not a proof that the
        // closed loop is resonant.
        //
        // So this returns the unresolved state, which is this file's stated
        // premise everywhere else: a quantity the arithmetic could not form has
        // produced no evidence, and an absence of evidence is neither a
        // certificate nor a refutation. See `riccati_accuracy` and
        // `resolved_magnitude` above for the same rule at the two other sites.
        // The closed-loop spectrum check in the discrete solver's verification
        // remains the SOLE owner of the stabilizing claim; this estimator makes
        // no claim about the spectrum at all.
        return {Scalar{0}, false};
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
/// P = U21 * U11^-1, symmetrize in place, and validate positive semi-definiteness.
/// Writes directly into P_out so that callers can avoid the ctrlpp::expected<Matrix>
/// return-by-value copy on the hot path.
///
/// The implementation prefers a back-substitution against U11 (solving
/// U11^T * P^T = U21^T) over an explicit inverse for numerical accuracy.
///
/// The positive semi-definiteness floor is the one psd_pivot_floor computes:
/// the floor factor times the state dimension times unit roundoff times the
/// largest absolute entry of P, negated. The operand is the largest entry and
/// NOT a row-sum norm -- the two differ by up to a factor of the dimension --
/// and the dimension factor is load-bearing rather than decorative. See that
/// function's docblock for the backward-error argument and for the measured
/// solution that missed the dimension-free floor while being positive
/// semi-definite to within one ulp; it is not restated here.
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
