#include "ctrlpp/control/dare.h"

#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <algorithm>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // Need 88 bytes: A(2x2=32) + B(2x1=16) + Q(2x2=32) + R(1x1=8)
    if(size < 88)
        return 0;

    double buf[11];
    std::memcpy(buf, data, 88);

    // Reject non-finite inputs early.
    for(int i = 0; i < 11; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    // Flush near-zero entries to exact zero. A raw magnitude far below the
    // clamp scale is kept as a tiny nonzero value it can build a uniformly-tiny
    // (and so falsely "well-scaled" relative to itself) controllability or
    // conditioning structure that a pure singular-value RATIO test does not
    // reject (two commensurately tiny singular values still divide out to a
    // benign-looking ratio), even though every entry is negligible against the
    // O(1)-O(2) problem scale this fuzzer explores. The floor is one part in a
    // hundred of the clamp bound, i.e. any entry more than two orders of
    // magnitude below the problem's own scale is treated as exact zero rather
    // than as a genuine small value.
    constexpr double b_clamp_bound = 2.0;
    const double zero_floor = b_clamp_bound / 1e2;
    auto clamp_entry = [zero_floor](double x) -> double
    {
        const double c = std::clamp(x, -2.0, 2.0);
        return std::abs(c) < zero_floor ? 0.0 : c;
    };

    // Clamp matrix entries to prevent intermediate overflow in symplectic construction
    Eigen::Matrix<double, 2, 2> A;
    A << clamp_entry(buf[0]), clamp_entry(buf[1]),
         clamp_entry(buf[2]), clamp_entry(buf[3]);

    // Reject a poorly-conditioned A via its own condition number. ctrlpp::dare
    // builds AinvT = A^-T as an intermediate step (see dare.h, Laub 1979 Eq. 7);
    // a near-singular A blows up AinvT's magnitude and propagates through the
    // whole symplectic construction. Same order-of-magnitude threshold and SVD
    // idiom as the controllability-matrix check below.
    // The smallest singular value is also checked against zero_floor directly
    // (not just the ratio to the largest): two commensurately tiny singular
    // values still divide out to a benign-looking condition number even though
    // the whole matrix is negligible against the problem's O(1)-O(2) scale.
    Eigen::JacobiSVD<Eigen::Matrix<double, 2, 2>> a_svd(A);
    const auto& a_sv = a_svd.singularValues();
    if(a_sv(1) < zero_floor || a_sv(0) / a_sv(1) > 10.0)
        return 0;

    Eigen::Matrix<double, 2, 1> B;
    B << clamp_entry(buf[4]), clamp_entry(buf[5]);

    Eigen::Matrix<double, 2, 2> Q_raw;
    Q_raw << clamp_entry(buf[6]), clamp_entry(buf[7]),
             clamp_entry(buf[8]), clamp_entry(buf[9]);

    // Reject a rank-deficient Q_raw: DARE's stabilizing solution is only
    // uniquely well-posed under detectability of (A, sqrt(Q)), and a full-rank
    // (invertible) Q_raw makes Q = Q_raw^T Q_raw strictly positive definite,
    // which trivially satisfies detectability for any A (an invertible output
    // map observes the full state directly). A rank-deficient Q_raw can make Q
    // singular in a direction that happens to align with an unstable mode of
    // A, giving a large-but-finite P whose absolute residual error grows with
    // its own magnitude; same rejection idiom as the controllability check
    // above.
    if(!Q_raw.colPivHouseholderQr().isInvertible())
        return 0;

    // Make Q positive semi-definite: Q = Q_raw^T * Q_raw
    Eigen::Matrix<double, 2, 2> Q = Q_raw.transpose() * Q_raw;

    // Make R positive definite: R = R_raw^2 + floor. The floor is derived from
    // B's own clamp bound (2.0) rather than eps: it keeps the worst-case
    // cheap-control amplification factor B^2/R (used in the tolerance below,
    // at B's clamp bound B^2 = 8) under a four-figure ceiling instead of
    // blowing up toward the reciprocal of machine epsilon, which is the
    // extreme near-singular regime where a direct Schur solve's absolute
    // (not relative) error necessarily grows without the oracle needing an
    // ever-larger, ad hoc safety margin to keep pace.
    double r_raw = std::clamp(buf[10], -2.0, 2.0);
    const double R_floor = (b_clamp_bound * b_clamp_bound) / 1e3;
    double R_val = r_raw * r_raw + R_floor;
    Eigen::Matrix<double, 1, 1> R;
    R << R_val;

    // Reject poorly-conditioned (A, B) pairs via the controllability matrix's
    // own condition number (largest / smallest singular value of [B, AB]),
    // rather than a binary invertibility test: the DARE stabilizing solution's
    // sensitivity to rounding grows with controllability conditioning (the
    // closed-loop gain K = S^-1 B'PA that drives the residual below scales
    // with it directly), so a moderately ill-conditioned pair -- not just an
    // exactly singular one -- can already blow up the absolute residual far
    // beyond what a fixed safety margin can absorb without becoming toothless.
    // The threshold is a single order of magnitude, comfortably above the
    // condition numbers observed on well-behaved random controllable pairs.
    Eigen::Matrix<double, 2, 2> ctrb;
    ctrb.col(0) = B;
    ctrb.col(1) = A * B;
    // As with A above, the smallest singular value is also checked directly
    // against zero_floor, not just its ratio to the largest.
    Eigen::JacobiSVD<Eigen::Matrix<double, 2, 2>> ctrb_svd(ctrb);
    const auto& ctrb_sv = ctrb_svd.singularValues();
    if(ctrb_sv(1) < zero_floor || ctrb_sv(0) / ctrb_sv(1) > 10.0)
        return 0;

    // Reject near-defective A: a small discriminant of the 2x2 characteristic
    // polynomial (trace^2 - 4*det) means A sits close to a non-diagonalizable
    // Jordan form. This is a well-known, independent source of ill-conditioning
    // for any Schur-based eigenstructure algorithm (clustered eigenvalues give
    // a poorly-conditioned eigenvector basis), unrelated to stabilizability;
    // the threshold is one percent of A's own squared Frobenius norm.
    const double trace_a = A.trace();
    const double det_a = A.determinant();
    const double discriminant = trace_a * trace_a - 4.0 * det_a;
    if(std::abs(discriminant) < A.squaredNorm() / 100.0)
        return 0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);

    // The library correctly declining on an ill-posed input (non-stabilizable,
    // singular, non-finite intermediate, etc.) is not a property violation.
    if(!result.has_value())
        return 0;

    const auto& P = result->P;

    // A finite-declared success must actually be finite.
    if(!P.allFinite())
        abort();

    // Live residual oracle: A'PA - P - A'PB (R + B'PB)^-1 B'PA + Q ~ 0.
    Eigen::Matrix<double, 2, 2> AtPA = A.transpose() * P * A;
    Eigen::Matrix<double, 1, 1> S = R + B.transpose() * P * B;
    Eigen::Matrix<double, 2, 2> cross = A.transpose() * P * B * S.inverse() * B.transpose() * P * A;
    Eigen::Matrix<double, 2, 2> resid = AtPA - P - cross + Q;

    // Backward-error tolerance: the residual is a sum of four n x n terms, each
    // formed from a chain of matrix products; the standard floating-point
    // matrix-multiplication backward-error bound accumulates rounding error
    // proportional to n per term, so the sum of the terms' own norms (the
    // "term_scale" below) is the honest base scale rather than just the raw
    // input norms -- this matters because AtPA and cross frequently sit at a
    // much larger common magnitude than their difference (the residual itself),
    // i.e. this is a catastrophic-cancellation regime, and a tolerance based on
    // the canceled result's own size would be far too tight. ctrlpp::dare's own
    // construction additionally forms G = B R^-1 B^T and AinvT = A^-T (see
    // dare.h, Laub 1979 Eq. 7) before the Schur decomposition of the symplectic
    // Z, so G's scale B^2/R is folded into the same term-magnitude sum. The
    // overall multiplicative margin below was calibrated empirically: run
    // against an extended fuzzing session (hundreds of thousands of random
    // controllable, well-scaled inputs) with no false-positive abort on the
    // library's correct default DARE path, following the same "eps-and-norm-
    // normalized residual vs. a generous integer-multiple threshold" style LAPACK's
    // own eigenvalue/Schur test suite (e.g. dget*, ddrvst) uses to absorb the
    // constant factors backward-error theory leaves unspecified.
    constexpr double lapack_style_margin = 20000.0;
    const double control_weight_scale = B.squaredNorm() / R(0, 0);
    const double term_scale = AtPA.norm() + P.norm() + cross.norm() + Q.norm() + control_weight_scale;
    const double n = static_cast<double>(P.rows());
    const double tol = lapack_style_margin
        * std::numeric_limits<double>::epsilon() * n * n * n * std::max(1.0, term_scale);

    if(resid.norm() > tol)
        abort();

    // P must be positive semi-definite: LDLT pivot-sign check against the same
    // floor the library's own extraction and postcondition use. The bound is
    // called rather than re-spelled here -- the two hand-copies of it that this
    // line and its fuzz_care twin used to carry stayed at `1 * eps * max|P_ij|`
    // when the library moved to the backward-error-carrying `N * eps *
    // max|P_ij|`, which made this oracle abort on a P that is positive
    // semi-definite to within one ulp.
    Eigen::LDLT<Eigen::Matrix<double, 2, 2>> ldlt(P);
    if(ldlt.info() != Eigen::Success
       || ldlt.vectorD().minCoeff() < ctrlpp::detail::psd_pivot_floor<double, 2>(P))
        abort();

    return 0;
}
