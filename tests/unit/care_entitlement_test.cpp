#include "fuzz/care_entitlement.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>

#include <bit>
#include <cmath>
#include <limits>
#include <optional>

namespace
{

using ctrlpp::fuzz::quad_reference_result;

struct pose
{
    double A[2][2];
    double B[2];
    double Q[2][2];
    double R;
};

// The fuzz target's budget: four times the ceil(log2(113)) Newton-Kleinman steps that double one
// correct bit up to binary128's 113-bit significand.
constexpr int reference_significand_bits = 113;
constexpr int reference_step_budget      = 4 * std::bit_width(unsigned{reference_significand_bits - 1});

constexpr double eps = std::numeric_limits<double>::epsilon();

// The accuracy-gate input of tests/fuzz/corpus/README.md as the target decodes it, with the
// target's Q = Q_raw' * Q_raw.
auto near_axis_pose() -> pose
{
    Eigen::Matrix2d Q_raw;
    Q_raw << -0x1.fbfbf31bfbf00p-4, -0x1.fbfbfbfbfbfbfp-4, 0.0, 0.0;
    const Eigen::Matrix2d Q = Q_raw.transpose() * Q_raw;
    return pose{{{0.0, 0.0}, {-2.0, -2.0}}, {-0x1.fbdbfbfbebfbfp-4, 0.0}, {{Q(0, 0), Q(0, 1)}, {Q(1, 0), Q(1, 1)}}, 0x1.004189374bc6ap+2};
}

// The three entry equations of A'P + PA - PBB'P + I = 0 give p12 = 1, p22 = sqrt(3) and p11 = p12 p22 at
// the stabilizing solution.
constexpr pose double_integrator{{{0.0, 1.0}, {0.0, 0.0}}, {0.0, 1.0}, {{1.0, 0.0}, {0.0, 1.0}}, 1.0};

auto refine(const pose &p, const double seed[2][2]) -> quad_reference_result
{
    return ctrlpp::fuzz::quad_refine_care(p.A, p.B, p.Q, p.R, seed, reference_step_budget);
}

auto refine_double_integrator() -> quad_reference_result
{
    const double root3      = std::sqrt(3.0);
    const double seed[2][2] = {{root3, 1.0}, {1.0, root3}};
    return refine(double_integrator, seed);
}

auto refine_near_axis(const pose &p) -> quad_reference_result
{
    // The reference's own answer from an identity seed, rounded to binary64; it only starts the iteration.
    const double seed[2][2] = {{0x1.f8115fc0b9661p-9, 0x1.f7ff7c8208fb3p-9}, {0x1.f7ff7c8208fb3p-9, 0x1.f7ff8e270cbffp-9}};
    return refine(p, seed);
}

auto condition_of(const pose &p, const quad_reference_result &reference) -> std::optional<double>
{
    return ctrlpp::fuzz::care_condition_number(p.A, p.B, p.Q, p.R, reference);
}

auto entitled(const pose &p, const quad_reference_result &reference) -> bool
{
    return ctrlpp::fuzz::care_forward_error_entitled(p.A, p.B, p.Q, p.R, reference);
}

auto entitled(const pose &p, const quad_reference_result &reference, double backward_error_factor) -> bool
{
    return ctrlpp::fuzz::care_forward_error_entitled(p.A, p.B, p.Q, p.R, reference, backward_error_factor);
}

}

TEST_CASE("The near-axis accuracy-gate pose is not entitled to a half-significand forward error", "[care][entitlement]")
{
    const pose p                          = near_axis_pose();
    const quad_reference_result reference = refine_near_axis(p);
    REQUIRE(reference.converged);

    const std::optional<double> condition = condition_of(p, reference);
    REQUIRE(condition.has_value());
    CHECK(*condition * eps > std::sqrt(eps));
    CHECK_FALSE(entitled(p, reference));
}

TEST_CASE("The double integrator is entitled to a half-significand forward error", "[care][entitlement]")
{
    const quad_reference_result reference = refine_double_integrator();
    REQUIRE(reference.converged);

    const std::optional<double> condition = condition_of(double_integrator, reference);
    REQUIRE(condition.has_value());
    CHECK(*condition * eps <= std::sqrt(eps));
    CHECK(entitled(double_integrator, reference));
}

TEST_CASE("An unconverged reference or a zero solution yields no condition number and no entitlement", "[care][entitlement]")
{
    SECTION("an unconverged reference")
    {
        quad_reference_result reference = refine_double_integrator();
        REQUIRE(reference.converged);
        reference.converged = false;

        CHECK_FALSE(condition_of(double_integrator, reference).has_value());
        CHECK_FALSE(entitled(double_integrator, reference));
    }

    SECTION("a converged zero solution")
    {
        constexpr pose unweighted{{{-1.0, 0.0}, {0.0, -1.0}}, {1.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}}, 1.0};
        const double seed[2][2]               = {{0.0, 0.0}, {0.0, 0.0}};
        const quad_reference_result reference = refine(unweighted, seed);
        REQUIRE(reference.converged);

        CHECK_FALSE(condition_of(unweighted, reference).has_value());
        CHECK_FALSE(entitled(unweighted, reference));
    }
}

TEST_CASE("Omitting the backward-error factor gives the verdict of a unit factor", "[care][entitlement]")
{
    const pose near_axis                             = near_axis_pose();
    const quad_reference_result near_axis_reference  = refine_near_axis(near_axis);
    const quad_reference_result integrator_reference = refine_double_integrator();
    REQUIRE(near_axis_reference.converged);
    REQUIRE(integrator_reference.converged);

    CHECK(entitled(near_axis, near_axis_reference) == entitled(near_axis, near_axis_reference, 1.0));
    CHECK(entitled(double_integrator, integrator_reference) == entitled(double_integrator, integrator_reference, 1.0));
}

TEST_CASE("A backward-error factor that carries the double integrator past the criterion withdraws its entitlement", "[care][entitlement]")
{
    const quad_reference_result reference = refine_double_integrator();
    REQUIRE(reference.converged);
    const std::optional<double> condition = condition_of(double_integrator, reference);
    REQUIRE(condition.has_value());

    const double factor = 2.0 / (*condition * std::sqrt(eps));
    CHECK_FALSE(entitled(double_integrator, reference, factor));
}
