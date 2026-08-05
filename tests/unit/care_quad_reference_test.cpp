#include "fuzz/care_quad_reference.h"

#include <catch2/catch_test_macros.hpp>

#include <cmath>

namespace
{

using ctrlpp::fuzz::quad;

// The double integrator, and WHY its solution is known in closed form. A reader
// must be able to check this without trusting either the reference below or the
// solver the reference exists to judge.
//
// The continuous equation is A'P + PA - P B R^-1 B' P + Q = 0. With
// A = [[0, 1], [0, 0]], B = [0, 1]', Q = I and R = 1, write P = [[p11, p12],
// [p12, p22]] and every term is available by inspection:
//
//   A'P + PA           = [[0, p11], [p11, 2*p12]]
//   P B R^-1 B' P      = [[p12^2, p12*p22], [p12*p22, p22^2]]
//
// so the equation is three scalar equations, one per distinct entry:
//
//   (1,1):  0    - p12^2     + 1 = 0   ->  p12^2 = 1
//   (1,2):  p11  - p12*p22       = 0   ->  p11   = p12 * p22
//   (2,2):  2*p12 - p22^2    + 1 = 0   ->  p22^2 = 2*p12 + 1
//
// The stabilizing solution is the positive semi-definite one, which fixes both
// signs: p12 = 1, p22 = sqrt(3) and p11 = p22. The other real symmetric
// solution, p22 = -sqrt(3), solves the same three equations and is used below as
// a seed OUTSIDE the stabilizing basin.
constexpr double pose_A[2][2] = {{0.0, 1.0}, {0.0, 0.0}};
constexpr double pose_B[2]    = {0.0, 1.0};
constexpr double pose_Q[2][2] = {{1.0, 0.0}, {0.0, 1.0}};
constexpr double pose_R       = 1.0;

// The step budget, counted rather than chosen. Newton-Kleinman doubles the
// number of correct significand bits per step inside its basin, so reaching a
// 113-bit significand from a seed carrying a single correct bit takes
// ceil(log2(113)) = 7 steps; the factor of four covers the approach from a
// stabilizing but inaccurate gain, before the quadratic rate takes hold. The
// budget bounds COST and not correctness -- under-sizing it produces
// abstentions, never wrong answers -- and the sweep below measures the counts it
// actually needs rather than assuming the margin.
constexpr int reference_significand_bits = 113;
constexpr int quadratic_steps_to_full_precision = []
{
    int steps = 0;
    for(int bits = 1; bits < reference_significand_bits; bits *= 2)
        ++steps;
    return steps;
}();
constexpr int pre_basin_approach_factor = 4;
constexpr int reference_step_budget     = pre_basin_approach_factor * quadratic_steps_to_full_precision;

auto refine_double_integrator(const double seed[2][2]) -> ctrlpp::fuzz::quad_reference_result
{
    return ctrlpp::fuzz::quad_refine_care(pose_A, pose_B, pose_Q, pose_R, seed, reference_step_budget);
}

/// @brief Round the reference's answer back to binary64.
///
/// The claim under test is that this rounding lands on the closed-form value bit
/// for bit, so the rounding is performed once, here, and compared exactly.
auto round_to_binary64(const quad P[2][2], double out[2][2]) -> void
{
    for(int i = 0; i < 2; ++i)
        for(int j = 0; j < 2; ++j)
            out[i][j] = static_cast<double>(P[i][j]);
}

/// @brief The closed loop a returned solution implies, on this pose.
///
/// Formed here from the entries rather than read off the reference's own
/// verdict, so that a guard which silently stopped firing would still be caught.
/// K = R^-1 B' P = [p12, p22] and A_cl = A - B K = [[0, 1], [-p12, -p22]].
auto closed_loop_of(const double P[2][2], double& trace, double& determinant) -> void
{
    const double gain[2] = {P[0][1] / pose_R, P[1][1] / pose_R};

    double closed_loop[2][2];
    for(int i = 0; i < 2; ++i)
        for(int j = 0; j < 2; ++j)
            closed_loop[i][j] = pose_A[i][j] - pose_B[i] * gain[j];

    trace       = closed_loop[0][0] + closed_loop[1][1];
    determinant = closed_loop[0][0] * closed_loop[1][1] - closed_loop[0][1] * closed_loop[1][0];
}

}

TEST_CASE("continuous quad reference reproduces a closed-form Riccati solution",
          "[care][reference]")
{
    const double root_three   = std::sqrt(3.0);
    const double seed[2][2]   = {{root_three, 1.0}, {1.0, root_three}};
    const auto   result       = refine_double_integrator(seed);

    REQUIRE(result.converged);

    double P[2][2];
    round_to_binary64(result.P, P);

    // The reference's own answer carries digits binary64 cannot hold: it sits
    // roughly 1e-16 off the binary64 grid, which is what makes the equalities
    // below a claim about rounding rather than a tautology.
    //
    // The CAPTURE is on the base-ten exponent because Catch2 stringifies a
    // double in FIXED notation, which would render this quantity as "0.0" and
    // tell a reader diagnosing a future failure the one thing that is not true
    // of it.
    const double grid_gap         = static_cast<double>(result.P[0][0] - quad{P[0][0]});
    const double grid_gap_decades = std::log10(std::abs(grid_gap));
    CAPTURE(result.steps, grid_gap_decades);

    // Bit equality, not an approximate comparison. The claim being made is
    // reproduction to the last bit of binary64, and a tolerance would quietly
    // weaken it into a different and much smaller claim. std::sqrt is correctly
    // rounded by IEEE 754, so the right-hand side is the binary64 value nearest
    // the exact root, which is what the reference must round onto.
    CHECK(P[0][0] == root_three);
    CHECK(P[0][1] == 1.0);
    CHECK(P[1][0] == 1.0);
    CHECK(P[1][1] == root_three);

    // The three properties the reference guards internally, re-derived here from
    // what it returned. Positive semi-definiteness of a real symmetric 2-by-2 is
    // all of its principal minors non-negative: the two diagonal entries and the
    // determinant.
    CHECK(P[0][1] == P[1][0]);
    CHECK(P[0][0] >= 0.0);
    CHECK(P[1][1] >= 0.0);
    CHECK(P[0][0] * P[1][1] - P[0][1] * P[1][0] >= 0.0);

    // Strict stability of the closed loop, as the Routh-Hurwitz condition on the
    // characteristic polynomial lambda^2 - trace*lambda + determinant: both roots
    // lie strictly left of the imaginary axis exactly when the trace is strictly
    // negative and the determinant strictly positive. No eigensolver and no
    // tolerance participate.
    double trace       = 0.0;
    double determinant = 0.0;
    closed_loop_of(P, trace, determinant);
    CHECK(trace < 0.0);
    CHECK(determinant > 0.0);
}

TEST_CASE("continuous quad reference converges from a perturbed seed inside its budget",
          "[care][reference]")
{
    const double root_three = std::sqrt(3.0);

    // Twelve decades of relative seed error, both signs. What is asserted is the
    // RELATION and not a count: strictly positive, strictly inside the budget,
    // and never rising as the seed improves -- which is the quadratic local rate
    // expressed as a step count. Pinning a count would assert the toolchain; the
    // observed counts belong in the record, not in a comparison.
    int previous_steps = reference_step_budget;

    for(int exponent = -1; exponent >= -12; --exponent)
    {
        for(int sign = 1; sign >= -1; sign -= 2)
        {
            const double relative = static_cast<double>(sign) * std::pow(10.0, static_cast<double>(exponent));
            const double factor   = 1.0 + relative;

            const double seed[2][2] = {{root_three * factor, factor},
                                       {factor, root_three * factor}};
            const auto   result     = refine_double_integrator(seed);

            CAPTURE(exponent, sign, result.steps, previous_steps);

            REQUIRE(result.converged);

            double P[2][2];
            round_to_binary64(result.P, P);

            CHECK(P[0][0] == root_three);
            CHECK(P[0][1] == 1.0);
            CHECK(P[1][0] == 1.0);
            CHECK(P[1][1] == root_three);

            CHECK(result.steps > 0);
            CHECK(result.steps < reference_step_budget);
            CHECK(result.steps <= previous_steps);

            previous_steps = result.steps;
        }
    }
}

TEST_CASE("continuous quad reference abstains rather than answering outside the basin",
          "[care][reference]")
{
    SECTION("a seed whose implied gain does not stabilize the pose")
    {
        // On this pose K = R^-1 B' P = [p12, p22], so A_cl = [[0, 1], [-p12, -p22]]
        // is stable only for p12 > 0 and p22 > 0. Seeded from the OTHER real
        // symmetric solution of the same equation, the iteration converges
        // cleanly -- it is already at a fixed point -- and returns a solution
        // that is negative definite. It is the definiteness guard, not the
        // convergence test, that withdraws the verdict here.
        //
        // That is the same guard that catches a dropped negation on the Lyapunov
        // right-hand side, which is the one transliteration error in this port
        // that does not announce itself: the wrong sign returns the negated
        // solution and converges just as cleanly as the right one.
        const double root_three = std::sqrt(3.0);
        const double seed[2][2] = {{-root_three, 1.0}, {1.0, -root_three}};
        const auto   result     = refine_double_integrator(seed);

        CAPTURE(result.steps);
        CHECK_FALSE(result.converged);
    }

    SECTION("a pose whose unstable mode is invisible to the weight")
    {
        // Zero weight and an unstable A. X = 0 solves the equation exactly, is
        // perfectly symmetric and is positive semi-definite, so the first two
        // guards pass; the closed loop it leaves is A itself, whose trace is +2.
        // Only the stability guard can withdraw this one, which is what makes it
        // a control on that guard rather than on the other two.
        const double A[2][2]    = {{1.0, 0.0}, {0.0, 1.0}};
        const double B[2]       = {0.0, 1.0};
        const double Q[2][2]    = {{0.0, 0.0}, {0.0, 0.0}};
        const double seed[2][2] = {{0.0, 0.0}, {0.0, 0.0}};

        const auto result = ctrlpp::fuzz::quad_refine_care(A, B, Q, 1.0, seed, reference_step_budget);

        CAPTURE(result.steps);
        CHECK_FALSE(result.converged);
    }

    SECTION("a closed loop on the imaginary axis")
    {
        // The continuous Lyapunov operator is singular whenever two closed-loop
        // eigenvalues sum to zero, which every purely imaginary conjugate pair
        // meets. The undamped oscillator with a zero seed leaves A_cl = A, whose
        // spectrum is {+i, -i}, so the operator is exactly singular and the
        // elimination finds a zero pivot. This is the failure mode the discrete
        // twin does not have: the Stein operator is singular only when a pair of
        // eigenvalues MULTIPLIES to one.
        const double A[2][2]    = {{0.0, 1.0}, {-1.0, 0.0}};
        const double B[2]       = {0.0, 1.0};
        const double Q[2][2]    = {{1.0, 0.0}, {0.0, 1.0}};
        const double seed[2][2] = {{0.0, 0.0}, {0.0, 0.0}};

        const auto result = ctrlpp::fuzz::quad_refine_care(A, B, Q, 1.0, seed, reference_step_budget);

        CAPTURE(result.steps);
        CHECK_FALSE(result.converged);
    }
}
