#include "ctrlpp/control/care.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>
#include <random>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <type_traits>

namespace
{

constexpr std::uint64_t sweep_seed = 0xc0a57ab1e5ULL;
constexpr int decade_count = 18;
constexpr int cases_per_decade = 64;
constexpr int cases_per_scale_direction = 32;

template <typename Scalar>
using matrix2 = Eigen::Matrix<Scalar, 2, 2>;

template <typename Scalar>
auto closed_loop_is_stable(const matrix2<Scalar>& A,
                           const matrix2<Scalar>& B,
                           const matrix2<Scalar>& R,
                           const matrix2<Scalar>& P) -> bool
{
    const matrix2<Scalar> closed_loop =
        (A - B * R.inverse() * B.transpose() * P).eval();
    Eigen::EigenSolver<matrix2<Scalar>> eigensystem(closed_loop, false);
    for(int index = 0; index < 2; ++index)
    {
        if(!(eigensystem.eigenvalues()(index).real() < Scalar{0}))
            return false;
    }
    return true;
}

// The margin BOTH gain comparisons in this file are held to.
//
// A relative gain error above the square root of epsilon is a wrong answer:
// for a radix-2 type that is exactly the loss of more than half the fractional
// significand bits, so it is derived from the type's radix rather than fitted,
// and it is the same criterion the discrete solver is judged by. It is a
// property of the answer rather than of the draw sequence.
//
// One margin, two comparisons, and they establish different things:
//
//  * against the same pose at unit weight scale -- a SCALE-INVARIANCE
//    CONSISTENCY CHECK. The continuous gain K = R^-1 B^T P does not move under
//    a common positive rescale of Q and R, so the two must agree. But the
//    reference is the same solver under the same method tag admitted by the
//    same acceptance rule, so an error common to both scales cannot be seen
//    this way, and a weight-scale defect IS that error.
//  * against the closed-form gain of the structurally simple family -- an
//    INDEPENDENT CRITERION, which owes the solver nothing.
//
// The two are compared against the same margin deliberately, so a difference
// between the columns is attributable to the reference rather than to the
// threshold.
const double wrong_gain_error =
    std::sqrt(std::numeric_limits<double>::epsilon());

// What the solver is REQUIRED to carry, and why it is this half of the sweep.
//
// A common rescale of Q and R by c leaves the closed-loop spectrum and the gain
// unchanged and multiplies the solution by c, so the residual test is
// scale-invariant. What is not scale-invariant is the extraction: the stable
// invariant subspace is spanned by [I; cP], so raising c drives the leading
// block U11 toward singularity and the extracted solution's error upward with
// it, while lowering c drives the subspace toward span[I; 0], which is the
// best-conditioned position it has. The two directions are therefore NOT
// symmetric, and only the downward one admits a guarantee the anchor can state
// without fitting a number.
//
// Downward it is stated and asserted as an equality: a rescale that does not
// increase the weights must not cost the solver its answer, over every decade
// swept. The one floor that direction ever hit -- the Householder tail whose
// square falls below the smallest normal value, which returned a zero solution
// as a success -- is now lifted by the equivariant rescale in the sign path, so
// the guarantee holds all the way down.
//
// Upward there is no such guarantee to state. Bisecting the accept/decline
// boundary for the identity family gives 1.0e+03 at n = 2, 1.3e+04 at n = 4 and
// 2.8e+04 at n = 8 -- neither order-one against the solver's counted rounding
// budget (ratios 9.3, 15.0 and 4.1) nor monotone in the dimension, and the
// two-dimensional comfortable family is not even monotone in the scale
// (declined at 1e3.0, accepted again at 1e3.6). A band asserted there would be
// a fitted constant, so acceptance is deliberately not asserted upward. What IS
// asserted upward is the property that matters: every answer the solver does
// accept is right, judged by the scale-invariant gain oracle above.
//
// This replaces an earlier equality over weights within a factor 1/epsilon of
// unit dynamics. That band's justification -- separability in double arithmetic
// -- was refuted by measurement: it forced the acceptance of 435 draws that the
// gain oracle shows to be wrong, so it was asserting the presence of wrong
// answers rather than the absence of over-rejection. Carrying the upper half of
// that band needs the weight equilibration the continuous path does not yet
// perform; until it does, the honest contract is a decline.

auto feedback_gain(const matrix2<double>& R,
                   const matrix2<double>& B,
                   const matrix2<double>& P) -> matrix2<double>
{
    return (R.inverse() * B.transpose() * P).eval();
}

// The relative distance between two gains, or NOTHING.
//
// The bare-scalar version of this helper divided by the reference magnitude
// unless that magnitude was zero, in which case it returned an UNNORMALIZED
// distance -- while every caller compared the return against a relative
// half-significand margin. A zero reference silently converted a relative test
// into an absolute one at whatever scale the draw happened to sit at, which is
// the degenerate comparison the resolved-magnitude primitive exists to
// eliminate.
//
// It now returns that primitive. A zero or non-finite reference yields the
// unresolved state, and `magnitude_within` refuses to read an unresolved
// operand as evidence of smallness, so the caller counts the pose as
// uncompared rather than as compared-and-passed.
auto relative_gain_error(const matrix2<double>& gain,
                         const matrix2<double>& reference)
    -> ctrlpp::detail::resolved_magnitude<double>
{
    const auto reference_scale = ctrlpp::detail::resolve_magnitude(reference);
    const auto difference = ctrlpp::detail::resolve_magnitude(
        (gain - reference).eval());
    if(!reference_scale.resolved || !difference.resolved
       || !(reference_scale.value > 0.0))
        return {0.0, false};
    return {difference.value / reference_scale.value, true};
}

// The gain the structurally simple family has EXACTLY, computed from the
// equation rather than from the solver.
//
// That family is A = -I, B = I, Q = sI, R = sI for a positive common weight
// multiple s. Substituting into A'P + PA - P B R^-1 B' P + Q = 0:
//
//     -P - P - P (1/s) P + sI = 0
//
// Every term is a multiple of the identity, so P = p I and the matrix equation
// collapses to one scalar quadratic:
//
//     p^2 + 2 s p - s^2 = 0   with p > 0
//
// whose positive root is p = s (sqrt(2) - 1). The gain follows:
//
//     K = R^-1 B' P = (1/s) * s (sqrt(2) - 1) I = (sqrt(2) - 1) I
//
// THE MULTIPLE CANCELS. The gain of this family is the same matrix at every
// weight scale, and it is a number this file can write down. Nothing in the
// derivation calls the solver, uses its acceptance rule, or reads its
// extraction code, which is exactly what the scale-invariance check below
// cannot say for itself.
auto structurally_simple_exact_gain() -> matrix2<double>
{
    return ((std::numbers::sqrt2 - 1.0) * matrix2<double>::Identity()).eval();
}

auto is_enumerated_care_error(ctrlpp::care_error error) -> bool
{
    switch(error)
    {
    case ctrlpp::care_error::non_lhp_stabilizable:
    case ctrlpp::care_error::non_finite_input:
    case ctrlpp::care_error::singular_r:
    case ctrlpp::care_error::singular_u11:
    case ctrlpp::care_error::non_psd_solution:
    case ctrlpp::care_error::schur_failed:
    case ctrlpp::care_error::sign_function_stagnated:
    case ctrlpp::care_error::unverified_solution:
        return true;
    }
    return false;
}

auto rotation(double angle) -> matrix2<double>
{
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    matrix2<double> result;
    result << cosine, -sine, sine, cosine;
    return result;
}

// One member of a swept population, built once and handed to every method tag.
//
// The draws are built ahead of the solving rather than interleaved with it so
// that the three public method tags are compared on literally the same
// matrices. Interleaving would make each tag's family a function of the order
// the generator was consumed in, and a difference between two tags would then
// be unattributable to either.
struct sweep_draw
{
    matrix2<double> A;
    matrix2<double> B;
    matrix2<double> Q;
    matrix2<double> R;
    // The same pose at unit weight scale. K = R^-1 B^T P does not move under a
    // common rescale of Q and R, so comparing the draw's gain against this
    // pose's gain is a SCALE-INVARIANCE CONSISTENCY CHECK.
    //
    // It is not an oracle for the acceptance rule, and it used to say it was.
    // The reference is produced by the same function, under the same method
    // tag, through the same extraction code, and is admitted by the same
    // acceptance rule as the answer it is being compared against. An error
    // COMMON TO BOTH SCALES is therefore invisible to it by construction -- and
    // a weight-scale defect is exactly that error class. What it does establish
    // is real and worth asserting: that the gain the solver reports does not
    // move when the weights are rescaled, which it is entitled to.
    //
    // The criterion that owes the solver nothing is the closed-form gain, and
    // it exists only for the structurally simple family. Its column is counted
    // separately below for that reason, rather than merged into this one.
    matrix2<double> unit_Q;
    matrix2<double> unit_R;
    bool            has_gain_oracle;
    bool            weights_reduced;
    // Whether this draw's gain is known in closed form, independently of the
    // solver. Only the structurally simple family is, and its value is
    // `structurally_simple_exact_gain()` at every weight scale.
    bool            has_closed_form_gain;
};

struct sweep_families
{
    std::vector<sweep_draw> near_axis;
    std::vector<sweep_draw> comfortable;
    std::vector<sweep_draw> structurally_simple;
};

auto build_sweep_families() -> sweep_families
{
    std::mt19937_64 generator(sweep_seed);
    std::uniform_real_distribution<double> mantissa(1.0, 10.0);
    std::uniform_real_distribution<double> alpha_distribution(0.2, 0.95);
    std::uniform_real_distribution<double> angle_distribution(
        -std::numbers::pi, std::numbers::pi);

    sweep_families families;

    for(int decade = 1; decade <= decade_count; ++decade)
    {
        const double decade_scale = std::pow(10.0, -decade);
        for(int index = 0; index < cases_per_decade; ++index)
        {
            const double delta = decade_scale * mantissa(generator);
            const double alpha = alpha_distribution(generator);
            const double beta = std::sqrt(1.0 - alpha * alpha);
            const double angle = angle_distribution(generator);
            const matrix2<double> transform = rotation(angle);

            matrix2<double> diagonal_a;
            diagonal_a << alpha * delta, 0.0, 0.0, -0.8;
            matrix2<double> diagonal_b;
            diagonal_b << beta * delta, 0.0, 0.0, 0.6;

            sweep_draw draw;
            draw.A = (transform * diagonal_a * transform.transpose()).eval();
            draw.B = (transform * diagonal_b).eval();
            draw.Q = matrix2<double>::Identity();
            draw.R = matrix2<double>::Identity();
            draw.unit_Q = draw.Q;
            draw.unit_R = draw.R;
            // The near-axis family carries no common weight rescale, so it has
            // no scale-invariant reference pose distinct from itself.
            draw.has_gain_oracle = false;
            draw.weights_reduced = false;
            draw.has_closed_form_gain = false;
            families.near_axis.push_back(draw);
        }

        for(int direction : {-1, 1})
        {
            for(int index = 0;
                index < cases_per_scale_direction;
                ++index)
            {
                const double exponent =
                    static_cast<double>(direction * decade);
                const double scale =
                    std::pow(10.0, exponent) * mantissa(generator);
                const matrix2<double> transform =
                    rotation(angle_distribution(generator));

                matrix2<double> base_a;
                base_a << -1.1, 0.3, -0.2, -1.4;
                matrix2<double> base_b;
                base_b << 0.8, 0.1, -0.15, 0.65;
                matrix2<double> base_q;
                base_q << 1.0, 0.2, 0.2, 1.7;
                matrix2<double> base_r;
                base_r << 1.3, 0.1, 0.1, 0.9;

                sweep_draw comfortable;
                comfortable.A =
                    (transform * base_a * transform.transpose()).eval();
                comfortable.B =
                    (transform * base_b * transform.transpose()).eval();
                comfortable.Q =
                    (scale * transform * base_q * transform.transpose()).eval();
                comfortable.R =
                    (scale * transform * base_r * transform.transpose()).eval();
                comfortable.unit_Q =
                    (transform * base_q * transform.transpose()).eval();
                comfortable.unit_R =
                    (transform * base_r * transform.transpose()).eval();
                comfortable.has_gain_oracle = true;
                comfortable.weights_reduced = scale < 1.0;
                comfortable.has_closed_form_gain = false;
                families.comfortable.push_back(comfortable);

                sweep_draw simple;
                simple.A = -matrix2<double>::Identity();
                simple.B = matrix2<double>::Identity();
                simple.Q = scale * matrix2<double>::Identity();
                simple.R = scale * matrix2<double>::Identity();
                simple.unit_Q = matrix2<double>::Identity();
                simple.unit_R = matrix2<double>::Identity();
                simple.has_gain_oracle = true;
                simple.weights_reduced = scale < 1.0;
                simple.has_closed_form_gain = true;
                families.structurally_simple.push_back(simple);
            }
        }
    }

    return families;
}

// What one method tag did to one population. Every field is a count, so no
// loop can pass vacuously: a census whose `drawn` is zero fails the same
// assertions a census full of wrong answers would.
struct population_census
{
    std::size_t drawn = 0;
    std::size_t accepted = 0;
    std::size_t accepted_unstable = 0;
    std::size_t accepted_with_wrong_gain = 0;
    // Disagreements against the CLOSED-FORM gain, which the solver had no part
    // in producing. Kept separate from the column above so the independent
    // criterion's verdict and the consistency check's verdict are visible side
    // by side rather than merged into one number that means neither.
    std::size_t closed_form_accepted = 0;
    std::size_t closed_form_compared = 0;
    std::size_t closed_form_wrong_gain = 0;
    double closed_form_worst_error = 0.0;
    // Comparisons that could not be made at all, because a reference magnitude
    // was zero or left the range. Counted rather than absorbed: a population
    // that quietly stopped being compared would otherwise report zero wrong
    // answers for the same reason an empty one would.
    std::size_t gain_comparison_unresolved = 0;
    std::size_t declined_unenumerated = 0;
    std::size_t accepted_separation_unavailable = 0;
    std::size_t accepted_reorder_incomplete = 0;
    std::size_t reduced_drawn = 0;
    std::size_t reduced_accepted = 0;
    std::size_t increased_drawn = 0;
    std::size_t increased_accepted = 0;
};

template <typename Method>
auto solve_population(const std::vector<sweep_draw>& draws, Method method_tag)
    -> population_census
{
    population_census census;
    for(std::size_t index = 0; index < draws.size(); ++index)
    {
        const sweep_draw& draw = draws[index];
        ++census.drawn;
        if(draw.has_gain_oracle)
        {
            if(draw.weights_reduced)
                ++census.reduced_drawn;
            else
                ++census.increased_drawn;
        }

        const auto result = ctrlpp::care<double, 2, 2>(
            draw.A, draw.B, draw.Q, draw.R, method_tag);
        CAPTURE(index);
        if(!result.has_value())
        {
            // Whatever a tag refuses, it refuses through the typed channel.
            if(!is_enumerated_care_error(result.error()))
                ++census.declined_unenumerated;
            CHECK(is_enumerated_care_error(result.error()));
            continue;
        }

        ++census.accepted;
        if(draw.has_gain_oracle)
        {
            if(draw.weights_reduced)
                ++census.reduced_accepted;
            else
                ++census.increased_accepted;
        }
        if(std::isnan(result->subspace_separation))
            ++census.accepted_separation_unavailable;
        if(!result->reorder_complete)
            ++census.accepted_reorder_incomplete;

        const bool stable =
            closed_loop_is_stable(draw.A, draw.B, draw.R, result->P);
        census.accepted_unstable += !stable;
        CHECK(stable);

        const matrix2<double> solved_gain =
            feedback_gain(draw.R, draw.B, result->P);

        // THE INDEPENDENT COLUMN, FIRST, because it is the one that can
        // contradict the solver. It is computed from the closed form derived
        // above and calls nothing in the library.
        //
        // Its margin is the SAME half-significand criterion the consistency
        // check uses, deliberately. The two columns exist to be read against
        // each other, and a column whose threshold is a counted-operation bound
        // would be answering a strictly tighter question -- whether the answer
        // is as good as the operation count entitles it to be -- so a
        // difference between the columns would not be attributable to the
        // reference. The worst error actually observed is recorded instead of
        // being bounded, so the tighter question has a measured answer without
        // a fitted constant being asserted.
        if(draw.has_closed_form_gain)
        {
            ++census.closed_form_accepted;
            const auto closed_form_error =
                relative_gain_error(solved_gain,
                                    structurally_simple_exact_gain());
            if(!closed_form_error.resolved)
            {
                ++census.gain_comparison_unresolved;
            }
            else
            {
                ++census.closed_form_compared;
                if(closed_form_error.value > census.closed_form_worst_error)
                    census.closed_form_worst_error = closed_form_error.value;
                CAPTURE(closed_form_error.value);
                census.closed_form_wrong_gain +=
                    closed_form_error.value > wrong_gain_error;
                CHECK(closed_form_error.value <= wrong_gain_error);
            }
        }

        if(!draw.has_gain_oracle)
            continue;

        const auto reference = ctrlpp::care<double, 2, 2>(
            draw.A, draw.B, draw.unit_Q, draw.unit_R, method_tag);
        REQUIRE(reference.has_value());
        const auto gain_error = relative_gain_error(
            solved_gain,
            feedback_gain(draw.unit_R, draw.B, reference->P));
        if(!gain_error.resolved)
        {
            // No comparison was made. Counted rather than skipped silently,
            // because an uncompared pose reports zero wrong answers for the
            // same reason an unsolved one does.
            ++census.gain_comparison_unresolved;
            continue;
        }
        CAPTURE(gain_error.value);
        census.accepted_with_wrong_gain += gain_error.value > wrong_gain_error;
        CHECK(gain_error.value <= wrong_gain_error);
    }
    return census;
}

// Magnitude at which the continuous path's stable-subspace factorization stops
// carrying the answer, written as its derivation.
//
// The path builds its stable subspace by a rank-revealing Householder QR of the
// projector (I - sign(H)) / 2. Every Householder step compares the SQUARED norm
// of its column tail against an absolute floor -- the scalar type's smallest
// normal value -- and discards the reflection when the tail falls to or below
// it. That tail carries the projector's stable/unstable coupling entry, which
// for the scalar family A = -a, B = Q = R = 1 is the solution itself, of size
// 1 / (2a). The square of the coupling therefore reaches the smallest normal
// value at
//
//     a* = 1 / (2 * sqrt(smallest normal value))
//
// and from there upward the whole band used to return SUCCESS carrying P = 0,
// whose Riccati residual is exactly 1 while the stabilizing answer is an
// ordinary normal double.
//
// The edge is computed from the scalar type's own limits because it IS that
// derivation. Its decimal value would pin the anchor to binary64 and would
// state nothing about the mechanism that produced it; the same expression
// locates the edge in binary32, 135 decades away.
template <typename Scalar>
auto subspace_collapse_magnitude() -> Scalar
{
    return Scalar{1}
           / (Scalar{2} * std::sqrt(std::numeric_limits<Scalar>::min()));
}

// The scalar family's Hamiltonian is [[-a, -1], [-1, a]], whose Frobenius sum
// of squares is 2a^2 + 2. It leaves the top of the range at a = sqrt(max / 2),
// so each anchor below states which regime it probes by this derivation rather
// than by quoting a magnitude.
template <typename Scalar>
auto hamiltonian_magnitude_ceiling() -> Scalar
{
    return std::sqrt(std::numeric_limits<Scalar>::max() / Scalar{2});
}

template <typename Scalar>
auto scalar_family_solve(Scalar magnitude)
    -> ctrlpp::expected<ctrlpp::care_result<Scalar, 1>, ctrlpp::care_error>
{
    Eigen::Matrix<Scalar, 1, 1> A;
    A << -magnitude;
    Eigen::Matrix<Scalar, 1, 1> B;
    B << Scalar{1};
    Eigen::Matrix<Scalar, 1, 1> Q;
    Q << Scalar{1};
    Eigen::Matrix<Scalar, 1, 1> R;
    R << Scalar{1};
    return ctrlpp::care<Scalar, 1, 1>(A, B, Q, R);
}

// P solves P^2 + 2aP - 1 = 0 with P > 0, hence P = 1 / (a + sqrt(a^2 + 1)).
// a^2 leaves the top of the range inside this band while a itself does not, so
// the root is taken in the scale-free form a * sqrt(1 + a^-2).
template <typename Scalar>
auto scalar_family_exact_solution(Scalar magnitude) -> Scalar
{
    const Scalar reciprocal = Scalar{1} / magnitude;
    const Scalar root =
        magnitude * std::sqrt(Scalar{1} + reciprocal * reciprocal);
    return Scalar{1} / (magnitude + root);
}

// The extracted solution passes through the rank-revealing QR of the 2n-by-2n
// projector and a triangular solve against U11. That count is the solver's own
// `care_extraction_rounding_ops`, quoted here rather than respelled: this bound
// and the solver's residual floor must be counted at the same basis size, and
// two independent spellings of one count is how they came to differ by a factor
// of eight. The reference root above costs three more operations. Every
// operation is counted whether or not it rounds.
template <typename Scalar>
auto scalar_family_accuracy_bound(Scalar exact) -> Scalar
{
    constexpr int extraction_rounding_ops =
        ctrlpp::detail::care_extraction_rounding_ops<1>;
    constexpr int reference_rounding_ops = 3;
    return Scalar{extraction_rounding_ops + reference_rounding_ops}
           * std::numeric_limits<Scalar>::epsilon() * exact;
}

// Real part of the eigenvalue of `A` nearest the imaginary axis. A continuous
// Riccati problem has a stabilizing solution only if every mode ON that axis is
// visible through the weighting, so locating the nearest mode is the first half
// of deciding whether an input is well posed.
auto nearest_axis_mode_abscissa(const matrix2<double>& A) -> double
{
    Eigen::EigenSolver<matrix2<double>> eigensystem(A, false);
    double nearest = std::numeric_limits<double>::infinity();
    for(int index = 0; index < 2; ++index)
    {
        const double abscissa =
            std::abs(eigensystem.eigenvalues()(index).real());
        if(abscissa < nearest)
            nearest = abscissa;
    }
    return nearest;
}

// The second half: how visible that mode is through the weighting, as a
// dimensionless number the solver had no part in producing.
//
// The quantity is the Rayleigh quotient v* Q v of the nearest-axis mode over the
// weighting's own magnitude. Eigen returns unit-norm eigenvectors, so the
// numerator carries no normalization of its own, and dividing by ||Q|| makes the
// result independent of how the caller scaled the weighting.
//
// **It is computed from the FACTOR, and passing the assembled Q would return
// zero.** For the input below the quotient is twenty-four decades under the
// entries of Q, so evaluating v' (L' L) v as a product against the assembled
// matrix cancels it away completely and reports exactly 0.0 -- which would then
// pass a "below the boundary" test for the wrong reason and fail a "positive"
// one that is true. Grouping it as ||L v||^2 does the subtraction in L v, where
// the result is around 1e-13 and perfectly representable, and only then squares.
// Q real symmetric gives v* Q v = ||L Re(v)||^2 + ||L Im(v)||^2, so no complex
// matrix product is needed either.
//
// Its meaning is a resolution threshold, not a well-posedness one. A mode on the
// axis whose visibility is positive but far below the arithmetic's ability to
// distinguish it from zero leaves a solution that EXISTS and is not determined:
// the solver must decline, and no accepted answer at that visibility could be
// checked against anything.
auto nearest_axis_mode_visibility(const matrix2<double>& A,
                                  const matrix2<double>& weight_factor)
    -> double
{
    Eigen::EigenSolver<matrix2<double>> eigensystem(A, true);
    double nearest = std::numeric_limits<double>::infinity();
    int nearest_index = 0;
    for(int index = 0; index < 2; ++index)
    {
        const double abscissa =
            std::abs(eigensystem.eigenvalues()(index).real());
        if(abscissa < nearest)
        {
            nearest = abscissa;
            nearest_index = index;
        }
    }

    const auto mode = eigensystem.eigenvectors().col(nearest_index);
    const Eigen::Vector2d real_part = mode.real();
    const Eigen::Vector2d imag_part = mode.imag();
    const double quotient = (weight_factor * real_part).squaredNorm()
                          + (weight_factor * imag_part).squaredNorm();

    // The denominator carries no cancellation -- it is a norm of nonnegative
    // contributions -- so it is taken from the assembled weighting directly.
    const double weight_magnitude =
        (weight_factor.transpose() * weight_factor).norm();
    return quotient / weight_magnitude;
}

}

TEST_CASE("CARE sign iteration declines an unresolved regenerated input",
          "[care][convergence][error]")
{
    matrix2<double> A;
    A << 0x0p+0, 0x0p+0,
        -0xf.dfdfdfdfdfdf8p-7, -0xf.dfdfdfdfdfdf8p-7;

    Eigen::Matrix<double, 2, 1> B;
    B << -0xf.dfdfdfdfdfdf8p-7,
         -0xf.dfdfdfdfdfdf8p-7;

    matrix2<double> raw_weight;
    raw_weight << -0xf.dfdfdfdfffdf8p-7, -0xf.dfdfdfdfdfdf8p-7,
                  -0xf.dfdfdfdfdf98p-7,  -0xf.dfdfdfdfdfdf8p-7;
    const matrix2<double> Q =
        (raw_weight.transpose() * raw_weight).eval();

    Eigen::Matrix<double, 1, 1> R;
    R << 0x8.3126e978d4fep-11;

    const auto result = ctrlpp::care<double, 2, 1>(A, B, Q, R);

    // WHY the refusal is required, stated without reference to the solver.
    //
    // A carries a mode exactly on the imaginary axis, and (A, B) is
    // controllable, so a stabilizing solution exists provided that mode is
    // visible through Q -- and it is, but at a relative visibility of 4.2e-25,
    // sixteen and a half decades below the point where binary64 can resolve it
    // and nine below epsilon itself. The solution
    // therefore exists and is not determined, which is a refusal.
    //
    // The boundary is sqrt(epsilon) rather than epsilon because the mode enters
    // the Hamiltonian through Q and leaves it through a squared quantity: a
    // visibility of v is distinguishable from zero in the iterate only while
    // v^2 stays above the roundoff of the sums that formed it. Sweeping the
    // visibility across thirty-one decades puts the accept/refuse transition at
    // 1e-8, against sqrt(epsilon) = 1.5e-8.
    const double axis_abscissa = nearest_axis_mode_abscissa(A);
    const double visibility = nearest_axis_mode_visibility(A, raw_weight);

    // The comparisons below are on the values themselves. The CAPTURE is on the
    // base-ten exponent because Catch2 stringifies a double in FIXED notation,
    // which renders every quantity in this test as "0.0" and would tell a reader
    // diagnosing a future failure the one thing that is not true of it.
    const double visibility_decades = std::log10(visibility);
    const double boundary_decades =
        std::log10(std::sqrt(std::numeric_limits<double>::epsilon()));
    CAPTURE(axis_abscissa, visibility_decades, boundary_decades);

    // That a mode sits exactly on the axis is a property of the literals, and it
    // is asserted against them rather than against an eigensolver: A is lower
    // triangular, so its spectrum IS its diagonal, and the leading diagonal
    // entry is exactly zero. No arithmetic and no threshold participate.
    CHECK(A(0, 1) == 0.0);
    CHECK(A(0, 0) == 0.0);

    // The COMPUTED abscissa is a different quantity and is not exactly zero on
    // every toolchain -- clang with FMA contraction returns 2.8e-18 where g++
    // returns exactly 0. A backward-stable eigensolver returns the exact
    // eigenvalues of A + E with ||E|| on the order of the operand count times
    // epsilon times ||A||, and this file counts a 2-by-2 at n = 2 the same way
    // the solver counts its own operations. Observed worst case is 28x inside
    // the bound; asserting equality with zero here asserted a toolchain.
    constexpr int spectral_rounding_ops = 2;
    CHECK(axis_abscissa <= double{spectral_rounding_ops}
                               * std::numeric_limits<double>::epsilon()
                               * A.norm());

    CHECK(visibility > 0.0);
    CHECK(visibility
          < std::sqrt(std::numeric_limits<double>::epsilon()));

    REQUIRE_FALSE(result.has_value());

    // WHICH refusal is not asserted, because it is not a property of the input.
    // Below the resolution boundary the enumerator is decided by instruction
    // selection: this input alone yields sign_function_stagnated under g++,
    // non_psd_solution under clang with contraction enabled, and
    // non_lhp_stabilizable under Apple clang on arm64, and a thirty-point
    // three-ulp neighborhood of it produces all three. Pinning one of them
    // asserted a toolchain.
    //
    // What the input DOES determine is the two diagnoses it can never carry: R
    // is nonsingular and every entry of every operand is finite. Either of those
    // would be the solver misreading a well-formed input as a domain violation.
    CHECK(result.error() != ctrlpp::care_error::singular_r);
    CHECK(result.error() != ctrlpp::care_error::non_finite_input);
    CHECK(is_enumerated_care_error(result.error()));

    const auto reference = ctrlpp::care<long double, 2, 1>(
        A.cast<long double>(),
        B.cast<long double>(),
        Q.cast<long double>(),
        R.cast<long double>());
    CHECK_FALSE(reference.has_value());
}

TEST_CASE("CARE sign iteration accepts only stable near-axis solutions",
          "[care][convergence][sweep]")
{
    const sweep_families families = build_sweep_families();

    // Acceptance is asserted per draw, against the draw's own weight scale,
    // rather than per decade against a count. A per-decade count is a property
    // of the generator's draw sequence, not of the solver: libstdc++, libc++
    // and MSVC give `std::uniform_real_distribution` different sequences from
    // the same seed, so an equality on the count fails on two of the three
    // while the solver behaves identically.
    const ctrlpp::detail::sign_function_care_method sign_tag;
    const population_census near_axis =
        solve_population(families.near_axis, sign_tag);
    const population_census comfortable =
        solve_population(families.comfortable, sign_tag);
    const population_census simple =
        solve_population(families.structurally_simple, sign_tag);

    CAPTURE(sweep_seed,
            near_axis.drawn,
            near_axis.accepted,
            near_axis.accepted_unstable);
    CHECK(near_axis.drawn
          == static_cast<std::size_t>(decade_count * cases_per_decade));
    CHECK(near_axis.accepted > 0);
    CHECK(near_axis.accepted_unstable == 0);

    // The sign path has no swap phase, so every answer it returns reports its
    // pivot-ratio diagnostic as unavailable and its reorder as complete.
    CHECK(near_axis.accepted_separation_unavailable == near_axis.accepted);
    CHECK(near_axis.accepted_reorder_incomplete == 0);
    CHECK(comfortable.accepted_separation_unavailable == comfortable.accepted);
    CHECK(simple.accepted_separation_unavailable == simple.accepted);

    // The over-rejection guard, on the half of the sweep where a guarantee can
    // be derived rather than fitted: a common rescale that does not increase the
    // weights must not cost the solver its answer, at any of the eighteen
    // decades swept downward. It is an equality because every such draw is
    // required to succeed.
    //
    // Upward the count is captured but not asserted, for the reason given at the
    // head of this file. The property asserted there instead is the one the
    // solver owes a caller: nothing it accepts is wrong. That check runs per
    // draw in both directions, against the draw's own unit-scale pose, and the
    // count below is its aggregate.
    const std::size_t reduced_drawn =
        comfortable.reduced_drawn + simple.reduced_drawn;
    const std::size_t reduced_accepted =
        comfortable.reduced_accepted + simple.reduced_accepted;
    const std::size_t increased_drawn =
        comfortable.increased_drawn + simple.increased_drawn;
    const std::size_t increased_accepted =
        comfortable.increased_accepted + simple.increased_accepted;
    const std::size_t accepted_with_wrong_gain =
        comfortable.accepted_with_wrong_gain + simple.accepted_with_wrong_gain;

    CAPTURE(sweep_seed,
            wrong_gain_error,
            reduced_drawn,
            reduced_accepted,
            increased_drawn,
            increased_accepted,
            accepted_with_wrong_gain);
    CHECK(reduced_drawn > 0);
    CHECK(increased_drawn > 0);
    CHECK(increased_accepted > 0);
    CHECK(reduced_accepted == reduced_drawn);
    CHECK(accepted_with_wrong_gain == 0);
    CHECK(simple.closed_form_wrong_gain == 0);
    CHECK(simple.closed_form_compared == simple.closed_form_accepted);
    CHECK(simple.closed_form_accepted > 0);
    CHECK(comfortable.gain_comparison_unresolved == 0);
    CHECK(simple.gain_comparison_unresolved == 0);
    CHECK(comfortable.declined_unenumerated == 0);
    CHECK(simple.declined_unenumerated == 0);
}

// The promise the public result contract makes is not qualified by method tag,
// so the evidence for it must not be either. The same three populations are
// swept under all three selectable tags, and every tag is held to the same two
// properties: nothing it accepts destabilizes the closed loop, and nothing it
// accepts disagrees with the scale-invariant gain oracle.
//
// The two Schur tags are the reason this case exists. Replayed against the
// revision before this one, they returned 150 and 156 accepted-unstable answers
// on the near-axis family, out of 814 and 839 accepted, while every member of
// that family is controllable and detectable.
TEMPLATE_TEST_CASE("CARE keeps the same promise under every method tag",
                   "[care][convergence][method]",
                   ctrlpp::detail::sign_function_care_method,
                   ctrlpp::detail::schur_care_method,
                   ctrlpp::detail::balanced_schur_care_method)
{
    const sweep_families families = build_sweep_families();
    const TestType method_tag;

    const population_census near_axis =
        solve_population(families.near_axis, method_tag);
    const population_census comfortable =
        solve_population(families.comfortable, method_tag);
    const population_census simple =
        solve_population(families.structurally_simple, method_tag);

    // The censuses are asserted nonempty first, so a tag whose loop never ran
    // cannot report zero wrong answers.
    CAPTURE(sweep_seed,
            near_axis.drawn,
            near_axis.accepted,
            near_axis.accepted_unstable,
            comfortable.drawn,
            comfortable.accepted,
            comfortable.accepted_unstable,
            comfortable.accepted_with_wrong_gain,
            simple.drawn,
            simple.accepted,
            simple.accepted_unstable,
            simple.accepted_with_wrong_gain,
            simple.closed_form_accepted,
            simple.closed_form_compared,
            simple.closed_form_wrong_gain,
            simple.closed_form_worst_error,
            near_axis.gain_comparison_unresolved,
            comfortable.gain_comparison_unresolved,
            simple.gain_comparison_unresolved);

    const std::size_t population_size =
        static_cast<std::size_t>(decade_count * cases_per_decade);
    CHECK(near_axis.drawn == population_size);
    CHECK(comfortable.drawn == population_size);
    CHECK(simple.drawn == population_size);
    CHECK(near_axis.accepted > 0);
    CHECK(comfortable.accepted > 0);
    CHECK(simple.accepted > 0);

    CHECK(near_axis.accepted_unstable == 0);
    CHECK(comfortable.accepted_unstable == 0);
    CHECK(simple.accepted_unstable == 0);

    CHECK(comfortable.accepted_with_wrong_gain == 0);
    CHECK(simple.accepted_with_wrong_gain == 0);

    // The independent column, asserted separately from the consistency column
    // above. Non-vacuity first: a census that compared nothing would otherwise
    // report zero disagreements for the same reason an empty one would.
    CHECK(simple.closed_form_accepted > 0);
    CHECK(simple.closed_form_compared == simple.closed_form_accepted);
    CHECK(simple.closed_form_wrong_gain == 0);
    CHECK(near_axis.closed_form_accepted == 0);
    CHECK(comfortable.closed_form_accepted == 0);

    // Every comparison that was owed was made. This is what keeps the two
    // zero-disagreement assertions above from being satisfiable by a
    // population that quietly stopped being compared.
    CHECK(near_axis.gain_comparison_unresolved == 0);
    CHECK(comfortable.gain_comparison_unresolved == 0);
    CHECK(simple.gain_comparison_unresolved == 0);

    CHECK(near_axis.declined_unenumerated == 0);
    CHECK(comfortable.declined_unenumerated == 0);
    CHECK(simple.declined_unenumerated == 0);

    // THE FULL-ACCEPTANCE CLAIM, PINNED.
    //
    // Two documents state in the present tense that the balanced tag answers
    // every draw of both common-scale populations, and one of them recommends
    // selecting that tag on exactly that basis. Nothing pinned it: the sweep
    // asserted only that acceptance exceeds zero, and the number had been
    // measured before an unconditional postcondition was added to that path.
    // It is re-measured and it holds, so it is asserted rather than written
    // down.
    //
    // THIS IS DELIBERATELY THE OVER-REJECTION GUARD WITH THE LEAST MARGIN IN
    // THE SUITE. It admits no declines at all on these two families, so it is
    // the first assertion a future tightening of the acceptance rule will
    // break. That is what it is for. A failure here is the guard firing and
    // says the tightening cost the balanced tag answers it used to give; it is
    // not a flaky test, and the two documents must move with it.
    //
    // THE GUARD IS GONE, because the reason for it is gone. It used to hold for
    // the balanced tag alone, whose diagonal similarity of the Hamiltonian was
    // the only equilibration the library had; the other two declined roughly
    // half of both populations. The solver now equilibrates the weightings
    // itself before the Hamiltonian is built -- which is the same similarity by
    // a different route -- so every tag answers every draw of both common-scale
    // populations, and the assertion is made of all three.
    CHECK(comfortable.accepted == comfortable.drawn);
    CHECK(simple.accepted == simple.drawn);
}

TEST_CASE("CARE carries the scalar family below the subspace-collapse edge",
          "[care][convergence][magnitude]")
{
    // One binary octave below the edge: every Householder tail the basis forms
    // is still normal, and this point is correct today. It is anchored so the
    // repair above it cannot be paid for by losing it.
    const double magnitude = subspace_collapse_magnitude<double>() / 2.0;
    const double exact = scalar_family_exact_solution(magnitude);
    const auto result = scalar_family_solve(magnitude);

    CAPTURE(std::log10(magnitude), std::log10(exact));
    REQUIRE(result.has_value());
    CHECK(result->P(0, 0) != 0.0);
    CHECK(std::abs(result->P(0, 0) - exact)
          <= scalar_family_accuracy_bound(exact));
    CHECK(magnitude < hamiltonian_magnitude_ceiling<double>());
}

TEST_CASE("CARE never returns a zero solution inside the finite-magnitude "
          "subspace-collapse band",
          "[care][convergence][magnitude]")
{
    // One binary octave above the edge. Every magnitude on the acceptance path
    // is still finite here -- the Hamiltonian's own Frobenius sum of squares
    // has not left the range -- so nothing about this point is an overflow.
    // The projector is idempotent and the extracted answer is zero, which is
    // exactly the pair the projector verdict must no longer be allowed to
    // certify.
    const double magnitude = subspace_collapse_magnitude<double>() * 2.0;
    const double exact = scalar_family_exact_solution(magnitude);
    const auto result = scalar_family_solve(magnitude);

    CAPTURE(std::log10(magnitude), std::log10(exact));
    CHECK(magnitude < hamiltonian_magnitude_ceiling<double>());

    // THIS CASE MEANS THE BAND MUST BE ANSWERED, and now says only that.
    //
    // It used to permit either outcome and then require an answer two lines
    // later, so the permissive branch could never be reached: a decline failed
    // at the requirement before the branch's assertion could mean anything.
    // Dead code in a test is worse than dead code elsewhere, because a reader
    // takes it as a statement about what the case tolerates. The equivariant
    // rescale reaches this magnitude, so an answer is owed here rather than a
    // refusal. This is the lower side of the recovered region's boundary; the
    // upper side is asserted in the case below, where BOTH outcomes are
    // genuinely admissible and the branch there is therefore kept.
    REQUIRE(result.has_value());
    CHECK(result->P(0, 0) != 0.0);
    CHECK(std::abs(result->P(0, 0) - exact)
          <= scalar_family_accuracy_bound(exact));
}

TEST_CASE("CARE never returns a zero solution where the Hamiltonian's own "
          "magnitude has left the range",
          "[care][convergence][magnitude]")
{
    // Five binary octaves above the edge, which places the Hamiltonian's
    // Frobenius sum of squares past the top of the range. This is the regime
    // the recorded counterexample sits in.
    const double magnitude = subspace_collapse_magnitude<double>() * 32.0;
    const double exact = scalar_family_exact_solution(magnitude);
    const auto result = scalar_family_solve(magnitude);

    CAPTURE(std::log10(magnitude), std::log10(exact));
    CHECK(magnitude > hamiltonian_magnitude_ceiling<double>());
    if(result.has_value())
    {
        CHECK(result->P(0, 0) != 0.0);
        CHECK(std::abs(result->P(0, 0) - exact)
              <= scalar_family_accuracy_bound(exact));
    }
    else
    {
        CHECK(is_enumerated_care_error(result.error()));
    }

    // The upper side of the recovered region's boundary. The Newton step forms
    // its own inverse, whose off-diagonal entry is of order 1/(a^2 + 1) and so
    // reaches exactly zero above one over the square root of the smallest
    // subnormal value; from there the accepted iterate carries exactly half the
    // true coupling and no rescale of the projector can restore it. The
    // solution is not recovered there and must therefore be refused, not
    // returned.
    const double halving_magnitude =
        1.0 / std::sqrt(std::numeric_limits<double>::denorm_min());
    const double far_magnitude = halving_magnitude * 2.0;
    const double far_exact = scalar_family_exact_solution(far_magnitude);
    const auto far_result = scalar_family_solve(far_magnitude);

    CAPTURE(std::log10(halving_magnitude),
            std::log10(far_magnitude),
            std::log10(far_exact));
    REQUIRE_FALSE(far_result.has_value());
    CHECK(is_enumerated_care_error(far_result.error()));
}

TEST_CASE("CARE sign iteration preserves its diagnostic contract",
          "[care][convergence][diagnostic]")
{
    const matrix2<float> A = -matrix2<float>::Identity();
    const matrix2<float> B = matrix2<float>::Identity();
    const matrix2<float> Q = matrix2<float>::Identity();
    const matrix2<float> R = matrix2<float>::Identity();

    const auto result = ctrlpp::care<float, 2, 2>(A, B, Q, R);

    REQUIRE(result.has_value());
    CHECK(closed_loop_is_stable(A, B, R, result->P));
    CHECK(std::isnan(result->subspace_separation));
    CHECK(result->reorder_complete);
}
