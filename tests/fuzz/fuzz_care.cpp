#include "ctrlpp/control/care.h"

#include "care_quad_reference.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <cstdio>
#include <limits>
#include <complex>
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

    // Clamp matrix entries to prevent intermediate overflow in Hamiltonian construction
    Eigen::Matrix<double, 2, 2> A;
    A << clamp_entry(buf[0]), clamp_entry(buf[1]),
         clamp_entry(buf[2]), clamp_entry(buf[3]);

    Eigen::Matrix<double, 2, 1> B;
    B << clamp_entry(buf[4]), clamp_entry(buf[5]);

    Eigen::Matrix<double, 2, 2> Q_raw;
    Q_raw << clamp_entry(buf[6]), clamp_entry(buf[7]),
             clamp_entry(buf[8]), clamp_entry(buf[9]);

    // ENTITLEMENT: the weight factor must make every unstable and marginal mode
    // of A visible. The stabilizing solution is uniquely well posed only under
    // detectability of the pair (A, weight factor), so detectability is what the
    // target must establish before it is entitled to demand an answer at all --
    // and once it is established, the binary128 policy iteration below converges
    // to the SAME stabilizing solution the library claims to have found, which
    // is what earns the forward-error framing of the verdict.
    //
    // The test is a rank test on the stacked pencil [A - lambda*I ; Q_raw] at
    // every eigenvalue lambda of A with non-negative real part. The stack has
    // full column rank at every such lambda exactly when no mode on or right of
    // the imaginary axis is unobservable through the weight.
    //
    // This REPLACES an invertibility test on the factor. Invertibility is
    // sufficient for detectability and far from necessary -- an invertible
    // output map observes the whole state -- so requiring it rejected every
    // rank-deficient weight, which is the ordinary regulator case where the
    // weight is an output map's Gram factor. That case was not explored at all
    // before this filter; it is now, deliberately.
    //
    // THE OPERAND IS THE RAW FACTOR AND NEVER THE ASSEMBLED WEIGHT. The raw
    // factor is a valid square-root factor of the assembled weight, so all the
    // information this test needs already lives in it, and testing the factor
    // means never forming the product and never paying the halving of available
    // digits that squaring costs. The same fact was met from the other side in
    // this tree: a quadratic form of an eigenvector evaluated against an
    // ASSEMBLED weight cancelled away completely and returned exactly 0.0, while
    // the identical quantity regrouped as a squared norm against the FACTOR did
    // the subtraction where it is representable and returned 1e-13.
    //
    // The eigenvalues are taken as complex and the stack is factorized in
    // complex arithmetic. Of a conjugate pair only the conjugate with
    // non-negative real part is tested: the other is its mirror and has the same
    // rank. Stacking a real system twice the size instead was considered and
    // rejected -- it doubles both dimensions and both operation counts below and
    // buys nothing.
    //
    // THE THRESHOLD IS THE SUM OF TWO COUNTED TERMS AND CARRIES NO FREE
    // COEFFICIENT.
    //
    //   (1) The factorization's own backward error. A column-pivoted Householder
    //       QR decides numerical rank at min(rows, cols) * eps relative to its
    //       largest pivot; that is the default the installed linear-algebra
    //       library applies, and its own source attributes the formula to Higham
    //       and notes it is the same one its LDLT already carries. The largest
    //       pivot of a column-pivoted factorization IS the largest column norm
    //       of the operand, so this term is a threshold relative to the stack's
    //       own norm, available without a second pass over it. Here
    //       min(rows, cols) is the column count, i.e. the state dimension.
    //   (2) The eigenvalue's own backward error. The shift is COMPUTED, not
    //       exact. A backward-stable eigensolver returns the exact eigenvalues
    //       of a nearby matrix whose perturbation is on the order of the operand
    //       count times epsilon times the norm (Golub & Van Loan, Matrix
    //       Computations, 4th ed., Sec. 7.5), and that perturbation enters the
    //       stack through the shifted block, whose multiplier is the identity
    //       and so has unit operator norm. This is the same counted form the
    //       continuous convergence anchor uses for its spectral margin, counted
    //       at the state dimension the same way.
    constexpr int state_dimension = 2;
    constexpr int stack_rank_ops = state_dimension;
    constexpr int eigenvalue_rounding_ops = state_dimension;

    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> a_eigensystem(A, false);
    if(a_eigensystem.info() != Eigen::Success)
        return 0;

    // `eigenvalues()` returns the vector BY VALUE, so the spectrum is named
    // before an entry of it is read. Reading through the unnamed call deduces an
    // expression holding a reference into a temporary that dies at the end of
    // its own statement; that failure is silent, optimization-dependent, and has
    // already turned continuous-integration legs of this tree red.
    const Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>>::EigenvalueType a_spectrum =
        a_eigensystem.eigenvalues();
    const double a_norm = A.norm();

    for(int index = 0; index < state_dimension; ++index)
    {
        const std::complex<double> mode = a_spectrum(index);
        if(mode.real() < 0.0)
            continue;

        Eigen::Matrix<std::complex<double>, 4, 2> pencil;
        pencil.topRows(2)    = A.cast<std::complex<double>>();
        pencil.bottomRows(2) = Q_raw.cast<std::complex<double>>();
        for(int diagonal = 0; diagonal < state_dimension; ++diagonal)
            pencil(diagonal, diagonal) -= mode;

        const Eigen::ColPivHouseholderQR<Eigen::Matrix<std::complex<double>, 4, 2>> pencil_qr(pencil);
        const double rank_threshold =
            double{stack_rank_ops} * std::numeric_limits<double>::epsilon() * pencil_qr.maxPivot()
            + double{eigenvalue_rounding_ops} * std::numeric_limits<double>::epsilon() * a_norm;

        // The rank is counted here rather than read off the factorization's own
        // rank(), which can only apply a threshold RELATIVE to its largest pivot
        // and therefore cannot carry the second term at all.
        const Eigen::Matrix<std::complex<double>, 4, 2>& pencil_factor = pencil_qr.matrixQR();
        int pencil_rank = 0;
        for(int diagonal = 0; diagonal < state_dimension; ++diagonal)
        {
            if(std::abs(pencil_factor(diagonal, diagonal)) > rank_threshold)
                ++pencil_rank;
        }

        if(pencil_rank != state_dimension)
            return 0;
    }

    // Make Q positive semi-definite: Q = Q_raw^T * Q_raw
    Eigen::Matrix<double, 2, 2> Q = Q_raw.transpose() * Q_raw;

    // Make R positive definite: R = R_raw^2 + floor. The floor is derived from
    // B's own clamp bound (2.0) rather than eps: it keeps the worst-case
    // cheap-control amplification factor B^2/R (used in the tolerance below,
    // at B's clamp bound B^2 = 8) under a four-figure ceiling instead of
    // blowing up toward the reciprocal of machine epsilon, which is the
    // extreme near-singular regime where a direct Schur/sign-function solve's
    // absolute (not relative) error necessarily grows without the oracle
    // needing an ever-larger, ad hoc safety margin to keep pace.
    double r_raw = std::clamp(buf[10], -2.0, 2.0);
    const double R_floor = (b_clamp_bound * b_clamp_bound) / 1e3;
    double R_val = r_raw * r_raw + R_floor;
    Eigen::Matrix<double, 1, 1> R;
    R << R_val;

    // Reject poorly-conditioned (A, B) pairs via the controllability matrix's
    // own condition number (largest / smallest singular value of [B, AB]),
    // rather than a binary invertibility test: the CARE stabilizing solution's
    // sensitivity to rounding grows with controllability conditioning (the
    // closed-loop gain K = R^-1 B'P that drives the residual below scales with
    // it directly), so a moderately ill-conditioned pair -- not just an
    // exactly singular one -- can already blow up the absolute residual far
    // beyond what a fixed safety margin can absorb without becoming toothless.
    // The threshold is a single order of magnitude, comfortably above the
    // condition numbers observed on well-behaved random controllable pairs.
    Eigen::Matrix<double, 2, 2> ctrb;
    ctrb.col(0) = B;
    ctrb.col(1) = A * B;
    // The smallest singular value is also checked against zero_floor directly
    // (not just the ratio to the largest): two commensurately tiny singular
    // values still divide out to a benign-looking condition number even though
    // the whole matrix is negligible against the problem's O(1)-O(2) scale.
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

    auto result = ctrlpp::care<double, 2, 1>(A, B, Q, R);

    // The library correctly declining on an ill-posed input (non-stabilizable,
    // singular, non-finite intermediate, etc.) is not a property violation.
    if(!result.has_value())
        return 0;

    const auto& P = result->P;

    // A finite-declared success must actually be finite.
    if(!P.allFinite())
        abort();

    // THE ACCURACY ORACLE. It does not consult the component it is judging.
    //
    // The answer must retain more than half of binary64's significand, measured
    // as its relative distance to an independently computed solution of the same
    // pose: Newton-Kleinman policy iteration at binary128, solving a continuous
    // Lyapunov equation per step and seeded from the returned solution. No
    // Hamiltonian matrix, no Schur decomposition, no Eigen, and above all no
    // call to the library's own acceptance machinery.
    //
    // That last exclusion is the point. The library requires a verdict from that
    // machinery BEFORE it returns, so any oracle built on it is downstream of a
    // decision the library has already made and cannot contradict it.
    //
    // The reference is seeded with the RETURNED SOLUTION, so the iteration
    // starts from the gain that solution implies and converges to the nearby
    // exact solution of the same pose. That is what makes the distance below a
    // forward error rather than a residual. It is also why the detectability
    // filter above had to land first: without it there is no guarantee the
    // iteration converges to the same stabilizing solution the library claims to
    // have found, and this framing would be unearned.
    //
    // THE STEP BUDGET IS DERIVED, AND UNDER-SIZING IT CANNOT CAUSE A FALSE
    // ABORT. Newton-Kleinman doubles the number of correct bits per step inside
    // its basin, so reaching binary128's 113-bit significand from a seed with a
    // single correct bit takes ceil(log2(113)) = 7 steps; the factor of four
    // covers the pre-basin approach from a stabilizing but inaccurate gain. A
    // budget that were too small would produce ABSTENTIONS, never aborts, so
    // this constant bounds cost rather than correctness. THE OBSERVED MAXIMUM ON
    // THE CONTINUOUS SIDE IS 9 STEPS, over the 276,588 poses that reached this
    // oracle in a campaign of 1,200,000 executions across eight pinned seeds and
    // four seedless runs, so the budget stands at three times the deepest
    // iteration any of them needed and not one abstention on that population was
    // budget-limited. That figure is a property of the population a campaign
    // draws rather than a bound on the iteration, and it was measured on one
    // architecture under one compiler.
    constexpr int reference_significand_bits = 113;
    constexpr int quadratic_steps_to_full_precision = []
    {
        int steps = 0;
        for(int bits = 1; bits < reference_significand_bits; bits *= 2)
            ++steps;
        return steps;
    }();
    constexpr int reference_step_budget = 4 * quadratic_steps_to_full_precision;

    double A_reference[2][2];
    double B_reference[2];
    double Q_reference[2][2];
    double P_reference[2][2];
    for(int i = 0; i < 2; ++i)
    {
        B_reference[i] = B(i, 0);
        for(int j = 0; j < 2; ++j)
        {
            A_reference[i][j] = A(i, j);
            Q_reference[i][j] = Q(i, j);
            P_reference[i][j] = P(i, j);
        }
    }

    const auto reference = ctrlpp::fuzz::quad_refine_care(A_reference, B_reference, Q_reference, R(0, 0), P_reference, reference_step_budget);

    // A reference that has not converged is not a reference. Its distance to the
    // answer measures nothing, so this pose yields no verdict. An absence of
    // evidence is not evidence of a defect, and it is emphatically not an abort.
    // The same disposition covers a distance that does not resolve because the
    // reference's own magnitude is zero: there is no scale to be relative to.
    ctrlpp::fuzz::quad distance_squared{};
    if(reference.converged && ctrlpp::fuzz::quad_relative_distance_squared(P_reference, reference.P, distance_squared))
    {
        // The half-significand criterion, derived from the radix in one
        // movement. Retaining more than half of a radix-two type's fractional
        // significand bits means a relative forward error below two to the minus
        // half that width; for binary64 that is 2^-26, which is exactly the
        // square root of the type's own epsilon, 2^-52. Squaring both sides
        // removes the square root -- one at binary128 would be a library call
        // and so a new dependency -- and leaves the comparison below. Nothing
        // here is fitted, and NO CONDITIONING FACTOR APPEARS ANYWHERE.
        //
        // THAT HAS A CONSEQUENCE THIS TARGET DOES NOT YET HANDLE, AND SAYING SO
        // HERE IS THE POINT. A criterion with no conditioning factor is not
        // neutral about conditioning; it is unconditionally TIGHT. The forward
        // error of any backward-stable algorithm at binary64 is about kappa*eps,
        // so a half-significand answer stops being attainable by ANY such
        // algorithm once kappa exceeds 1/sqrt(eps), about 6.7e7. Past that point
        // this comparison aborts on an answer sitting at the arithmetic's own
        // floor. An infinite tolerance would be an oracle that cannot fail; an
        // unconditionally tight one is an oracle that fails on correct behavior,
        // and the two are the same mistake pointed in opposite directions.
        //
        // A campaign found such a pose rather than one being built to make the
        // point. A = [0 0 ; -2 -2], whose marginal mode at zero has eigenvector
        // [1, -1]/sqrt(2) and is visible through the weight factor at 3.3e-8,
        // giving a closed-loop abscissa of -2.05e-9 and a Lyapunov separation
        // whose reciprocal is 2.44e8 -- well past the 6.7e7 above. The returned
        // answer's relative forward error there is 3.5e-8 against a
        // backward-stable prediction of kappa*eps = 5.4e-8, so the answer is AT
        // the bound, and this comparison aborts at 2.33 times the criterion. The
        // reference converged in five steps, its own answer is positive
        // semi-definite, and both other gates pass, so none of that is a
        // reference artifact. Walking one weight entry across 41 units in the
        // last place at that same visibility, the solver declines on 37 rungs and
        // answers on 4, and all four fail here: at this conditioning the
        // acceptance decision is selected by rounding rather than by what is
        // attainable.
        //
        // WHAT IS MISSING IS AN ENTITLEMENT TEST, NOT A LOOSER CRITERION. The
        // filters above decide entitlement on the OPEN-loop state matrix and on
        // the pair (A, B); neither sees the closed loop's distance to the
        // imaginary axis, which is what governs this solution's conditioning. A
        // target that demands a half-significand answer has to establish first
        // that it is asking for one that exists in the arithmetic it asks in.
        // Choosing that test's form, and its constant, on the single pose above
        // would be fitting to a population of one, which is why it is recorded
        // here and not written here.
        const ctrlpp::fuzz::quad margin_squared{std::numeric_limits<double>::epsilon()};
        if(distance_squared > margin_squared)
            abort();
    }

    // THE RESIDUAL IS REPORTED AND NEVER COMPARED, AND THAT IS DELIBERATE.
    //
    // It is computed by a different arithmetic route than the library uses -- an
    // explicit inverse of R rather than the solve the library performs -- but on
    // the continuous side that buys nothing an oracle may act on. The discrete
    // twin can bound its own residual because it has two arithmetic routes into
    // ONE estimator and their agreement is a real property. Here there is only
    // one estimator of continuous accuracy, and it is precisely the one this
    // oracle must not consult. Any bound placed on this residual would therefore
    // be a second tolerance whose sole justification is that it is loose, which
    // is the exact defect this target carried and the reason its oracle was
    // rewritten. THE VALUE ESTABLISHES NOTHING ABOUT ACCURACY. The accuracy
    // criterion is the independent forward error above; this quantity exists so
    // a campaign can report its distribution beside that verdict.
    Eigen::Matrix<double, 2, 2> AtP = A.transpose() * P;
    Eigen::Matrix<double, 1, 1> Rinv = R.inverse();
    Eigen::Matrix<double, 2, 2> cross = P * B * Rinv * B.transpose() * P;
    Eigen::Matrix<double, 2, 2> resid = AtP + AtP.transpose() - cross + Q;

    // Reporting is off unless a campaign asks for it. This target processes
    // millions of inputs, so an unconditional line per input would be a cost
    // rather than a report; the lookup is done once. The reference's step count
    // rides along, because separating a budget-limited abstention from genuine
    // ill-conditioning needs the count and not just the verdict.
    static const bool report_residual = std::getenv("CTRLPP_FUZZ_REPORT_CARE_RESIDUAL") != nullptr;
    if(report_residual)
    {
        std::fprintf(stderr, "care residual %.17g reference_steps %d reference_converged %d\n",
                     resid.norm(), reference.steps, reference.converged ? 1 : 0);
    }

    // P must be positive semi-definite: LDLT pivot-sign check against the same
    // floor the library's own extraction uses, called rather than re-spelled so
    // the two cannot drift apart again (see the fuzz_dare twin).
    //
    // THIS CHECK IS WHY THIS TARGET IS LISTED AS KNOWN RED, AND IT IS NEITHER OF
    // THE TWO GATES ABOVE. What it refuses is characterized here rather than left
    // to be recovered from a recorded artifact, because the characterization is
    // what says a longer campaign cannot settle it.
    //
    //  * On a returned solution that is rank one to within rounding, this pivot
    //    is a Schur complement of two nearly equal quantities, so it takes only
    //    INTEGER multiples of one unit in the last place of the largest entry of
    //    P. Measured to five significant figures on both recorded refusals.
    //  * The floor is a RELATIVE quantity, floor_factor * N * eps * max|P_ij|.
    //    Its width in those same integer units runs from exactly two at the
    //    bottom of a binade to just under four at the top, because eps times a
    //    value is a fixed fraction of it while a unit in the last place is a
    //    step. So the check admits a pivot of one or two units and refuses at
    //    three -- but only while the largest entry sits in the lower half of its
    //    binade. The identical three-unit pivot on the identical matrix shape is
    //    admitted higher up: measured on the floor primitive alone, refused at
    //    1.05 and 1.25 and admitted at 1.50, 1.75 and 1.99. Both recorded
    //    refusals sit within five percent of the bottom, where it is tightest.
    //  * The refusal is therefore a BAND and not a knife edge. Walking one input
    //    field in units in the last place around a refusing pose refuses on
    //    thirteen consecutive rungs in one band and five in another, while a
    //    census of 298,112 poses drawn by a campaign of 1,200,000 executions
    //    produced none at all. The region is one a random draw essentially never
    //    enters and a plateau once entered, which is why refusals arrive in
    //    clusters or not at all, and why one clean campaign settles nothing.
    //  * THE ANSWERS INSIDE THAT BAND ARE ACCURATE. Across every swept rung the
    //    forward error squared sits between 5.4e-33 and 3.4e-28, eleven to
    //    sixteen decades inside the accuracy criterion above, and on six rungs
    //    the accuracy oracle ISSUED its verdict, the verdict PASSED, and this
    //    check refused the same solution anyway.
    //  * Six build configurations reproduce the whole band byte for byte --
    //    three optimization levels, two compilers, both linear-algebra patch
    //    releases and the sanitizer configuration this tree builds -- so it is a
    //    property of the arithmetic rather than of instruction selection. Fused
    //    multiply-add contraction removes every refusal, which is the mechanism
    //    confirming itself: the three-unit pivot exists because three separate
    //    roundings accumulate, and fusing them collapses it.
    //
    // THE FLOOR WAS NOT WIDENED, and its exposed factor is left at its default.
    // A coefficient chosen after observing pivots at three units is a constant
    // fitted to a measured population whatever derivation is written beside it.
    // The open question is whether the floor should carry a second term for the
    // rounding already present in P as the solver returned it, since today it
    // counts only the factorization's own. That is a question about a library
    // primitive both solvers' acceptance machinery calls rather than about this
    // target, which is why it is not answered here.
    Eigen::LDLT<Eigen::Matrix<double, 2, 2>> ldlt(P);
    if(ldlt.info() != Eigen::Success
       || ldlt.vectorD().minCoeff() < ctrlpp::detail::psd_pivot_floor<double, 2>(P))
        abort();

    // THE CLOSED-LOOP STABILITY CHECK, ONE-SIDED AND WITH A STATED RESOLUTION.
    //
    // A stabilizing solution must place every closed-loop eigenvalue strictly in
    // the open left half-plane. This is the one accuracy-independent property
    // that still means something on the draws where the reference yields no
    // verdict, which is why it is kept as its own check rather than folded into
    // the forward error above.
    //
    // The margin is counted, not chosen: a backward-stable eigensolver returns
    // the exact eigenvalues of a nearby matrix whose perturbation is on the
    // order of the operand count times epsilon times that matrix's norm (Golub &
    // Van Loan, Matrix Computations, 4th ed., Sec. 7.5), and this file counts a
    // two-by-two the same way the solver counts its own operations. A zero
    // margin was rejected: a computed abscissa is not exactly zero on every
    // toolchain, and asserting against zero asserts a toolchain rather than a
    // property.
    //
    // THE CHECK IS STRICTLY ONE-SIDED. It aborts only on a RESOLVED refutation,
    // a real part above the POSITIVE margin. Inside the band between the
    // negative and the positive margin it yields NO VERDICT. The near-defective
    // filter above tests the OPEN-LOOP state matrix, not this closed loop, so on
    // a near-defective closed loop the eigenvalue condition number amplifies the
    // backward error and the counted bound understates the uncertainty; a
    // two-sided assertion there would abort on a resolution this check does not
    // have. A margin exists to give the check a stated resolution, and reading
    // it honestly means declining inside it.
    //
    // The library's own quasi-triangular spectral predicate is deliberately not
    // called here, for the same reason its acceptance machinery is not: the
    // solver already required a verdict from it before returning.
    constexpr int spectral_rounding_ops = state_dimension;

    // Every operand is named before a member or a column of it is read. The
    // unnamed form deduces a block holding a reference into a returned temporary
    // that dies at the end of the statement, which is silent,
    // optimization-dependent and has already turned five continuous-integration
    // legs of this tree red.
    const Eigen::Matrix<double, 1, 2> gain = R.inverse() * B.transpose() * P;
    const Eigen::Matrix<double, 2, 2> closed_loop = A - B * gain;

    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> closed_loop_eigensystem(closed_loop, false);
    if(closed_loop_eigensystem.info() == Eigen::Success)
    {
        const Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>>::EigenvalueType closed_loop_spectrum =
            closed_loop_eigensystem.eigenvalues();
        const double stability_margin = double{spectral_rounding_ops}
                                        * std::numeric_limits<double>::epsilon()
                                        * closed_loop.norm();

        for(int index = 0; index < state_dimension; ++index)
        {
            if(closed_loop_spectrum(index).real() > stability_margin)
                abort();
        }
    }

    return 0;
}
