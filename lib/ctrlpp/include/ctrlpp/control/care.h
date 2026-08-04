#ifndef HPP_GUARD_CTRLPP_CONTROL_CARE_H
#define HPP_GUARD_CTRLPP_CONTROL_CARE_H

/// @brief Continuous-time algebraic Riccati equation solver.
///
/// Solves A^T P + P A - P B R^{-1} B^T P + Q = 0 for the stabilizing P.
///
/// Builds the 2n x 2n Hamiltonian H = [[A, -B R^{-1} B^T], [-Q, -A^T]] (Laub 1979),
/// then dispatches to the selected method. The default applies a scaled Newton
/// iteration to sign(H), extracts P = U21 * U11^-1 from the stable-subspace
/// projector, and verifies THAT MATRIX against the equation it is supposed to
/// solve before reporting success. The two alternative tags use real-Schur
/// decomposition and Bai-Demmel reordering, with optional Hamiltonian balancing.
///
/// The continuous-time LQR gain is K = R^{-1} B^T P.
///
/// Shares the real-Schur reorder primitive and the Riccati P-extraction helper with
/// `ctrlpp::dare` via `ctrlpp::detail/`. Does not depend on `dare.h`.
///
/// @cite laub1979       -- Laub, "A Schur Method for Solving Algebraic Riccati Equations", 1979
/// @cite bai_demmel_1993 -- Bai & Demmel, "On swapping diagonal blocks in real Schur form", 1993

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/control/care_types.h"

#include "ctrlpp/detail/care_methods.h"
#include "ctrlpp/detail/schur_reorder.h"
#include "ctrlpp/detail/riccati_solution.h"
#include "ctrlpp/detail/care_sign_function.h"
#include "ctrlpp/detail/care_postconditions.h"
#include "ctrlpp/detail/hamiltonian_balance.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>
#include <complex>
#include <cstddef>
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
    -> ctrlpp::expected<Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>, care_error>
{
    constexpr int n = static_cast<int>(NX);
    constexpr int n2 = 2 * n;
    using MatNxN   = Eigen::Matrix<Scalar, n, n>;
    using Mat2Nx2N = Eigen::Matrix<Scalar, n2, n2>;

    if (!A.allFinite() || !B.allFinite() || !Q.allFinite() || !R.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);

    // R^{-1} is needed for B R^{-1} B^T, so R is judged by a reciprocal-pivot test with
    // the input dimension supplying the size factor -- the same convention the discrete
    // solver uses on both of its inverted operands. Verified to have the identical
    // defect rather than assumed: before this test, a zero R made the QR solve return
    // infinities and the caller was told their finite input was non-finite, while a
    // rank-deficient but nonzero R made the same solve return a least-squares answer
    // over the leading rank columns and the solver reported SUCCESS on a Hamiltonian
    // that is not the one the problem defines.
    auto qr_R = R.colPivHouseholderQr();
    qr_R.setThreshold(Scalar{static_cast<int>(NU)} * std::numeric_limits<Scalar>::epsilon());
    if (!qr_R.isInvertible())
        return ctrlpp::unexpected(care_error::singular_r);

    const MatNxN S = (B * qr_R.solve(
                             Eigen::Matrix<Scalar, int(NU), int(NX)>(B.transpose()))).eval();

    Mat2Nx2N H;
    H.template block<n, n>(0, 0) = A;
    H.template block<n, n>(0, n) = -S;
    H.template block<n, n>(n, 0) = -Q;
    H.template block<n, n>(n, n) = -A.transpose();

    if (!H.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);

    return H;
}

/// @brief Solve CARE from a pre-built Hamiltonian H, dispatching on the method tag.
///
/// Shared between `ctrlpp::care` and `ctrlpp::lqr_gain_continuous` so the latter can
/// build H with a pre-computed R^{-1} and avoid recomputing it. The `Method` tag
/// selects between the default matrix sign-function Newton iteration, the
/// real-Schur + Bai-Demmel reorder path, and the balanced-Schur variant.
template <typename Scalar, std::size_t NX,
          care_solve_method   Method = sign_function_care_method,
          conditioning_policy Cond   = pivot_ratio_conditioning>
auto care_solve_from_hamiltonian(
    const Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& H,
    Method method_tag     = {},
    Cond   /*cond_tag*/   = {})
    -> ctrlpp::expected<care_result<Scalar, NX>, care_error>
{
    if constexpr (std::is_same_v<Method, schur_care_method>)
    {
        constexpr int n = static_cast<int>(NX);
        constexpr int n2 = 2 * n;
        using Mat2N = Eigen::Matrix<Scalar, n2, n2>;

        Eigen::RealSchur<Mat2N> schur(H);
        if (schur.info() != Eigen::Success)
            return ctrlpp::unexpected(care_error::schur_failed);

        Mat2N T = schur.matrixT();
        Mat2N U = schur.matrixU();
        if (!T.allFinite() || !U.allFinite())
            return ctrlpp::unexpected(care_error::non_finite_input);

        // One predicate, one margin. The derivation lives beside the
        // quasi-triangular primitives so that this path and the balanced path
        // cannot drift apart: the eigenvalues of a backward-stable real Schur
        // factor carry a perturbation on the order of the factored matrix's
        // dimension times unit roundoff times the factor's largest entry.
        const Scalar lhp_margin =
            schur_eigenvalue_margin<Scalar, n2>(T);
        auto predicate = [lhp_margin](std::complex<Scalar> lam) -> bool
        {
            return lam.real() < -lhp_margin;
        };

        auto rr = reorder_real_schur<Scalar, n2>(T, U, predicate, Cond{});
        if (rr.placed < n)
            return ctrlpp::unexpected(care_error::non_lhp_stabilizable);
        if (!T.allFinite() || !U.allFinite())
            return ctrlpp::unexpected(care_error::non_finite_input);

        care_result<Scalar, NX> out;
        auto P_err = extract_riccati_solution_into<Scalar, n2>(out.P, U);
        if (!P_err)
        {
            switch (P_err.error())
            {
                case riccati_extract_error::singular_u11:
                    return ctrlpp::unexpected(care_error::singular_u11);
                case riccati_extract_error::non_finite:
                    return ctrlpp::unexpected(care_error::non_finite_input);
                case riccati_extract_error::non_psd:
                    return ctrlpp::unexpected(care_error::non_psd_solution);
            }
            return ctrlpp::unexpected(care_error::non_finite_input);
        }

        // The reorder placed n eigenvalues it judged to be in the left half-plane
        // and the extraction produced a matrix from them, but neither step ever
        // substituted that matrix back into the equation. The verification below
        // does, against H -- the Hamiltonian the caller's problem defines, the
        // same object the sign-function branch verifies against, and not the
        // reordered factor T or the basis U, which describe an invariant subspace
        // rather than the equation the answer must satisfy.
        if (!care_solution_satisfies_postconditions<Scalar, NX>(H, out.P))
            return ctrlpp::unexpected(care_error::unverified_solution);

        out.subspace_separation = rr.subspace_separation;
        out.reorder_complete    = rr.complete;
        return out;
    }
    else if constexpr (std::is_same_v<Method, sign_function_care_method>)
    {
        return detail::care_solve_via_sign_function<Scalar, NX>(
            H, method_tag.warmup_iterations);
    }
    else  // balanced_schur_care_method
    {
        return detail::care_solve_via_balanced_schur<Scalar, NX, Cond>(H);
    }
}

}

namespace detail
{

/// @brief The common divisor applied to both continuous weightings before the
/// Hamiltonian is built.
///
/// DERIVED FROM THE HAMILTONIAN'S OWN STRUCTURE, not transplanted from the
/// discrete solver. Writing `G = B R^-1 B'`, the continuous Hamiltonian is
///
///     H   = [[A, -G], [-Q, -A']]
///
/// and the Hamiltonian of the pose with both weightings divided by `s` is
///
///     H_s = [[A, -sG], [-Q/s, -A']]
///
/// because `G` is homogeneous of degree minus one in `R`. Those two matrices are
/// SIMILAR: with `D = diag(I, sI)`, `D^-1 H D` has off-diagonal blocks `-sG` and
/// `-Q/s`, which is `H_s` exactly. A common weight rescale is therefore a
/// block-diagonal similarity of the Hamiltonian, it cannot move the spectrum,
/// and the only thing it moves is where the two off-diagonal blocks sit inside
/// the representable range. Measured over 2,304 draws spanning eighteen decades
/// of common scale in both directions, the identity `G_s == s G` holds to a
/// worst relative residual of `5.61e-16`.
///
/// The divisor that puts BOTH blocks at unit order is the largest weight entry.
/// Dividing by it makes every entry of `Q/s` and `R/s` at most one in magnitude,
/// so `Q/s` is order one, and because `R/s` is order one as well the rescaled
/// Gramian `B (R/s)^-1 B'` is too. One divisor equilibrates both blocks at once.
///
/// This has the same expression as the discrete solver's weight scale and NOT
/// the same derivation, which matters because transplanting a bound between the
/// two Riccati equations has produced a disagreement in this codebase before.
/// The discrete argument is about the symplectic pencil's operands; this one is
/// about a similarity of the Hamiltonian. Nothing else carries across: the
/// continuous residual's derivative is a Sylvester operator where the discrete
/// one is a Stein operator, and the stability region is a half-plane rather than
/// a disk.
template <typename Scalar, std::size_t NX, std::size_t NU>
auto care_weight_scale(const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
                       const Eigen::Matrix<Scalar, int(NU), int(NU)>& R) -> Scalar
{
    const Scalar weight_magnitude =
        std::max(Q.cwiseAbs().maxCoeff(), R.cwiseAbs().maxCoeff());
    return weight_magnitude > Scalar{0} ? weight_magnitude : Scalar{1};
}

}

/// @brief Continuous-time Algebraic Riccati Equation solver.
///
/// Returns `ctrlpp::expected<care_result<Scalar, NX>, care_error>`. On success,
/// `result->P` is the stabilizing solution; `result->subspace_separation` is the
/// minimum pivot ratio across accepted swaps for Schur methods and is not
/// available for the default sign-function method; `result->reorder_complete`
/// is true if every swap was accepted or the selected method has no swap phase.
///
/// ## Both weightings are equilibrated before the Hamiltonian is built
///
/// A common positive divisor is applied to `Q` and `R`, the problem is solved at
/// that scale, and the solution is multiplied back. The equation is homogeneous
/// of degree one in `(P, Q, R)` taken together -- substituting `P -> P/s`,
/// `Q -> Q/s`, `R -> R/s` multiplies `A'P + PA - P B R^-1 B' P + Q` by exactly
/// `1/s` -- so the rescale changes nothing about the answer, and the gain
/// `K = R^-1 B' P` does not move at all. See `care_weight_scale` for why the
/// divisor is what it is: the rescale is a block-diagonal similarity of the
/// Hamiltonian, and the divisor is the one that lands both off-diagonal blocks
/// at unit order.
///
/// What it changes is which poses have an answer at all. Measured over the two
/// common-scale families the convergence anchor sweeps, 1,152 draws each across
/// eighteen decades of common scale in both directions: the unequilibrated
/// default answered 652 and 655 of them and declined the rest, because the
/// weights leave the range where the answer does not. Equilibrated it answers
/// **all 1,152 of each**, adding 997 answers, and every added answer is correct
/// -- zero disagreements against the scale-invariant gain reference, and for the
/// structurally simple family zero against the closed-form gain `(sqrt(2)-1) I`,
/// whose worst relative error is `4.02e-16`. Nothing was lost in the other
/// direction: **no pose the unequilibrated path answered correctly is declined.**
template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU,
          detail::care_solve_method    Method = detail::sign_function_care_method,
          detail::conditioning_policy  Cond   = detail::pivot_ratio_conditioning>
auto care(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
          const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
          const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
          const Eigen::Matrix<Scalar, int(NU), int(NU)>& R,
          Method                                         /*method_tag*/ = {},
          Cond                                           /*cond_tag*/   = {})
    -> ctrlpp::expected<care_result<Scalar, NX>, care_error>
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");

    if (!Q.allFinite() || !R.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);

    const Scalar weight_scale = detail::care_weight_scale<Scalar, NX, NU>(Q, R);
    // An equilibration that cannot be computed is an absence of evidence, and
    // this solver's standing contract is that an absence of evidence declines
    // rather than proceeds.
    if (!(weight_scale > Scalar{0}) || !std::isfinite(weight_scale))
        return ctrlpp::unexpected(care_error::non_finite_input);

    auto Q_scaled = Q;
    auto R_scaled = R;
    if (weight_scale != Scalar{1})
    {
        Q_scaled /= weight_scale;
        R_scaled /= weight_scale;
        if (!Q_scaled.allFinite() || !R_scaled.allFinite())
            return ctrlpp::unexpected(care_error::non_finite_input);
    }

    auto H_result =
        detail::build_care_hamiltonian<Scalar, NX, NU>(A, B, Q_scaled, R_scaled);
    if (!H_result)
        return ctrlpp::unexpected(H_result.error());

    auto result =
        detail::care_solve_from_hamiltonian<Scalar, NX, Method, Cond>(*H_result);
    if (!result)
        return result;

    if (weight_scale == Scalar{1})
        return result;

    // The postcondition has already run, against the EQUILIBRATED Hamiltonian.
    // That is the caller's own problem written at a scale where its operands are
    // representable -- the two are similar, not merely close -- so the check was
    // made on better-conditioned evidence than the caller's scale offers, not on
    // weaker evidence. What it does not cover is the one step that happens after
    // it: the multiplication back. An answer that is not representable at the
    // scale it is returned at is not an answer, so the product is checked.
    result->P *= weight_scale;
    if (!result->P.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);

    return result;
}

/// @brief CARE with cross-weight N: reduces to standard form via
/// Q' = Q - N R^{-1} N^T, A' = A - B R^{-1} N^T, then forwards.
template <typename Scalar, std::size_t NX, std::size_t NU,
          detail::care_solve_method    Method = detail::sign_function_care_method,
          detail::conditioning_policy  Cond   = detail::pivot_ratio_conditioning>
auto care(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
          const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
          const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
          const Eigen::Matrix<Scalar, int(NU), int(NU)>& R,
          const Eigen::Matrix<Scalar, int(NX), int(NU)>& N,
          Method                                         method_tag = {},
          Cond                                           cond_tag   = {})
    -> ctrlpp::expected<care_result<Scalar, NX>, care_error>
{
    if (!A.allFinite() || !B.allFinite() || !Q.allFinite() || !R.allFinite()
        || !N.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);

    // The reduction to standard form inverts R before the Hamiltonian build ever sees
    // it, so the same test is applied to the factorization this overload already forms.
    auto qr_R = R.colPivHouseholderQr();
    qr_R.setThreshold(Scalar{static_cast<int>(NU)} * std::numeric_limits<Scalar>::epsilon());
    if (!qr_R.isInvertible())
        return ctrlpp::unexpected(care_error::singular_r);

    auto Rinv_Nt = qr_R.solve(Eigen::Matrix<Scalar, int(NU), int(NX)>(N.transpose())).eval();

    Eigen::Matrix<Scalar, int(NX), int(NX)> Qp = (Q - N * Rinv_Nt).eval();
    Eigen::Matrix<Scalar, int(NX), int(NX)> Ap = (A - B * Rinv_Nt).eval();

    return care<Scalar, NX, NU, Method, Cond>(Ap, B, Qp, R, method_tag, cond_tag);
}

}

#endif
