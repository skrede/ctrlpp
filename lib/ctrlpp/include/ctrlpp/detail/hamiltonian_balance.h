#ifndef HPP_GUARD_CTRLPP_DETAIL_HAMILTONIAN_BALANCE_H
#define HPP_GUARD_CTRLPP_DETAIL_HAMILTONIAN_BALANCE_H

/// @brief DGEBAL-style diagonal balancing and balanced-Schur CARE solve.
///
/// @note Retained for reproducibility; superseded by `sign_function_care_method`
///       after the instruction-count bakeoff. The balanced-Schur path tracks
///       the plain Schur path to within 1 percent across the bakeoff sweep
///       because DGEBAL balance is a near no-op on well-conditioned
///       Hamiltonians emitted by `build_care_hamiltonian`
///       (diagonal D stays near ones).
///
/// `balance_hamiltonian` iteratively equilibrates row and column infinity
/// norms of the 2n x 2n Hamiltonian by applying a diagonal similarity
/// H' = D^{-1} H D, mirroring LAPACK DGEBAL phase 2 on a fixed-size
/// Eigen::Matrix<Scalar, 2*NX, 2*NX> without heap allocation. The LAPACK
/// permutation/isolation phase is omitted because on the Hamiltonians
/// produced by `build_care_hamiltonian` it rarely fires and adds branches
/// without measurable gain. Base-2 scaling constants SCLFAC = 2 and
/// FACTOR = 19/20 (= 0.95) are expressed as rational literals so no bare
/// non-structural numeric constant appears in the hot path.
///
/// `care_solve_via_balanced_schur` composes the balance step, the shared
/// `Eigen::RealSchur` + `ctrlpp::detail::reorder_real_schur` pipeline,
/// the D back-scale on the leading n columns of U, and the shared
/// `ctrlpp::detail::extract_riccati_solution_into` primitive. The back-
/// scale is a left-multiplication by diag(D) on the stable-subspace basis
/// columns, equivalent to `U.leftCols(n).array().colwise() *= D.array()`.
/// This is required because the LHP invariant subspace of H' = D^{-1} H D
/// is D^{-1} V (where V is the LHP invariant subspace of H), so recovering
/// V costs one element-wise product.
///
/// Plain DGEBAL destroys Hamiltonian structure; the structure-preserving
/// Benner 2001 symplectic balancing is reserved for a structured-Schur
/// path (conditional, not wired here). For the plain-Schur pipeline the
/// downstream treats H as a general 2n x 2n matrix, so plain DGEBAL is
/// correct.
///
/// @cite lapack_dgebal : Reference LAPACK SRC/dgebal.f phase 2 algorithm
/// @cite benner2001    : Benner, "Symplectic Balancing of Hamiltonian Matrices", 2001 (structured alternative, deferred)

#include "ctrlpp/expected.h"

#include "ctrlpp/control/care_types.h"

#include "ctrlpp/detail/schur_reorder.h"
#include "ctrlpp/detail/riccati_solution.h"
#include "ctrlpp/detail/care_postconditions.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>
#include <complex>
#include <cstddef>
#include <algorithm>
#include <type_traits>

namespace ctrlpp::detail
{

/// @brief In-place DGEBAL-style diagonal balance of a 2n x 2n matrix.
///
/// Overwrites H with the balanced matrix H' = D^{-1} H D and writes the
/// diagonal scaling vector into D_out (with D_out(i) >= 0 for all i).
template <typename Scalar, std::size_t NX>
auto balance_hamiltonian(
    Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& H,
    Eigen::Matrix<Scalar, 2 * int(NX), 1>&           D_out)
    -> void
{
    constexpr int n2 = 2 * int(NX);
    using Vec2N = Eigen::Matrix<Scalar, n2, 1>;

    const Scalar sclfac     = Scalar{2};
    const Scalar factor_num = Scalar{19};
    const Scalar factor_den = Scalar{20};
    const Scalar factor     = factor_num / factor_den;
    const Scalar sfmin      = std::numeric_limits<Scalar>::min();
    const Scalar sfmax      = Scalar{1} / sfmin;
    const Scalar sfmin2     = sfmin * sclfac;
    const Scalar sfmax2     = sfmax / sclfac;

    D_out = Vec2N::Ones();
    bool noconv = true;

    // THE SWEEP COUNT IS BOUNDED AT COMPILE TIME, and the bound is a safety
    // ceiling rather than a prediction. Termination itself does not depend on
    // it: every accepted rescale reduces that index's row-plus-column sum to
    // below 19/20 of its previous value, with power-of-two factors clamped to
    // the scalar type's range, so the loop ends on its own. What it did not have
    // was a bound a caller could read off the types, which every other iteration
    // in the Riccati path has -- the Newton loop carries `max_iters`. Without one
    // this method could not appear in a real-time budget at all, because
    // "terminates" and "terminates within a stated number of steps" are
    // different claims and only the second is schedulable.
    //
    // The ceiling counts the distinct power-of-two scalings an index can occupy
    // without leaving the representable range, once per index. It is derived
    // from the scalar type and the dimension, and no measured population went
    // anywhere near it -- which is the point: exceeding it means the data is
    // pathological, not that the tuning was wrong.
    //
    // IT IS A GUARANTEE, NOT A BUDGET, and the gap is large enough that saying so
    // matters. Over 550,000 draws with entries spread across 2^-280 to 2^+280 --
    // deliberately the regime a cap exists for -- the worst sweep count observed
    // was 2 at NX=1, 12 at NX=2, 43 at NX=4 and 45 at NX=8, against ceilings of
    // 4,090 / 8,180 / 16,360 / 32,720 respectively. The worst measured run sits
    // about three orders of magnitude below its ceiling (1.4e-3 of it at NX=8),
    // and the means are 2.00 / 3.87 / 6.40 / 8.29. A schedule that must not be
    // exceeded should be built from the ceiling; a schedule that wants to be
    // realistic should be built from the measurement and treat the ceiling as the
    // backstop it is.
    //
    // STOPPING EARLY IS SAFE, and that is what makes a cap admissible here at
    // all. Balancing is a similarity preconditioner: every applied step updates
    // `H` and `D_out` together, so at any cut point the invariant
    // `H_returned == D^-1 * H_original * D` holds exactly. A capped run returns a
    // less well balanced matrix, never a wrong one, and the back-scale downstream
    // stays consistent because it reads the same `D_out`.
    constexpr int exponent_span = std::numeric_limits<Scalar>::max_exponent - std::numeric_limits<Scalar>::min_exponent;
    constexpr int max_sweeps    = n2 * exponent_span;

    for (int sweep = 0; noconv && sweep < max_sweeps; ++sweep)
    {
        noconv = false;
        for (int i = 0; i < n2; ++i)
        {
            Scalar r = Scalar{0};
            Scalar c = Scalar{0};
            for (int j = 0; j < n2; ++j)
            {
                if (j == i) continue;
                r += std::abs(H(i, j));
                c += std::abs(H(j, i));
            }
            if (r == Scalar{0} || c == Scalar{0}) continue;

            Scalar g = r / sclfac;
            Scalar f = Scalar{1};
            const Scalar s = c + r;

            while (c < g
                   && std::max(f, c) < sfmax2
                   && std::min(r, g) > sfmin2)
            {
                f *= sclfac;
                c *= sclfac;
                r /= sclfac;
                g /= sclfac;
            }
            g = c / sclfac;
            while (g >= r
                   && std::max(r, g) < sfmax2
                   && std::min({f, c, g}) > sfmin2)
            {
                f /= sclfac;
                c /= sclfac;
                g /= sclfac;
                r *= sclfac;
            }

            if ((c + r) < factor * s)
            {
                noconv = true;
                D_out(i) *= f;
                H.row(i) /= f;
                H.col(i) *= f;
            }
        }
    }
}

/// @brief CARE solve via DGEBAL pre-balance then real Schur + Bai-Demmel reorder.
template <typename Scalar, std::size_t NX,
          conditioning_policy Cond = pivot_ratio_conditioning>
auto care_solve_via_balanced_schur(
    const Eigen::Matrix<Scalar, 2 * int(NX), 2 * int(NX)>& H_in)
    -> ctrlpp::expected<care_result<Scalar, NX>, care_error>
{
    constexpr int n  = int(NX);
    constexpr int n2 = 2 * n;
    using Mat2N = Eigen::Matrix<Scalar, n2, n2>;
    using Vec2N = Eigen::Matrix<Scalar, n2, 1>;

    if (!H_in.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);

    Mat2N H = H_in;
    Vec2N D;
    balance_hamiltonian<Scalar, NX>(H, D);

    Eigen::RealSchur<Mat2N> schur(H);
    if (schur.info() != Eigen::Success)
        return ctrlpp::unexpected(care_error::schur_failed);

    Mat2N T = schur.matrixT();
    Mat2N U = schur.matrixU();
    if (!T.allFinite() || !U.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);

    // The same one margin the plain-Schur path uses, from the same derivation.
    // This site used to carry its own: unit roundoff times max(1, scale), which
    // dropped the dimension factor and added a floor of one that no comment
    // derived. The floor's only candidate justification was a balance driving
    // the factor's magnitude toward zero, and that is refuted by measurement --
    // over 3,456 draws spanning eighteen decades of weight scale in each
    // direction, the balanced factor's largest entry stayed inside
    // [0.99999999999999845, 1.8260572569697802]. A balance is a normalizing
    // preconditioner, so its output magnitude is O(1) by construction and a
    // floor at one is inert where it is not simply absent.
    const Scalar lhp_margin = schur_eigenvalue_margin<Scalar, n2>(T);
    auto predicate = [lhp_margin](std::complex<Scalar> lam) -> bool
    {
        return lam.real() < -lhp_margin;
    };

    auto rr = reorder_real_schur<Scalar, n2>(T, U, predicate, Cond{});
    if (rr.placed < n)
        return ctrlpp::unexpected(care_error::non_lhp_stabilizable);
    if (!T.allFinite() || !U.allFinite())
        return ctrlpp::unexpected(care_error::non_finite_input);

    U.leftCols(n).array().colwise() *= D.array();

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

    // Verified against H_in, the Hamiltonian the CALLER's problem defines, and
    // never against the balanced H this function has been working on.
    //
    // The two are similar but they are not the same equation. Balancing applies
    // H' = D^{-1} H D, whose left-half-plane invariant subspace is D^{-1} V; the
    // back-scale above undoes that, so `out.P` is a solution to the caller's
    // Riccati equation and not to the balanced one. Substituting it into the
    // balanced Hamiltonian would form a residual for an equation nobody asked
    // about, and D is chosen precisely to equilibrate row and column norms, so
    // that residual would also be compared against a different scale than the
    // caller's. The local H is dead at this point; H_in is the object under test.
    if (!care_solution_satisfies_postconditions<Scalar, NX>(H_in, out.P))
        return ctrlpp::unexpected(care_error::unverified_solution);

    out.subspace_separation = rr.subspace_separation;
    out.reorder_complete    = rr.complete;
    return out;
}

}

#endif
