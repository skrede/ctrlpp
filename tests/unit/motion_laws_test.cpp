#include "ctrlpp/traj/motion_laws.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <numbers>

using namespace ctrlpp;
using Catch::Matchers::WithinAbs;

// --- Cubic law ---

TEST_CASE("cubic_law at tau=0", "[traj][motion_laws]")
{
    auto pt = cubic_law(0.0);
    CHECK_THAT(pt.q, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(6.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(-12.0, 1e-12));
}

TEST_CASE("cubic_law at tau=1", "[traj][motion_laws]")
{
    auto pt = cubic_law(1.0);
    CHECK_THAT(pt.q, WithinAbs(1.0, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(-6.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(-12.0, 1e-12));
}

TEST_CASE("cubic_law at tau=0.5", "[traj][motion_laws]")
{
    auto pt = cubic_law(0.5);
    CHECK_THAT(pt.q, WithinAbs(0.5, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(1.5, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(-12.0, 1e-12));
}

// --- Quintic law ---

TEST_CASE("quintic_law at tau=0", "[traj][motion_laws]")
{
    auto pt = quintic_law(0.0);
    CHECK_THAT(pt.q, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(60.0, 1e-12));
}

TEST_CASE("quintic_law at tau=1", "[traj][motion_laws]")
{
    auto pt = quintic_law(1.0);
    CHECK_THAT(pt.q, WithinAbs(1.0, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(0.0, 1e-12));
    // dddq(1) = 60 - 360 + 360 = 60 (symmetric with tau=0 value)
    CHECK_THAT(pt.dddq, WithinAbs(60.0, 1e-12));
}

TEST_CASE("quintic_law at tau=0.5 (symmetric midpoint)", "[traj][motion_laws]")
{
    auto pt = quintic_law(0.5);
    CHECK_THAT(pt.q, WithinAbs(0.5, 1e-12));
}

// --- Septic law ---

TEST_CASE("septic_law at tau=0", "[traj][motion_laws]")
{
    auto pt = septic_law(0.0);
    CHECK_THAT(pt.q, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(0.0, 1e-12));
}

TEST_CASE("septic_law at tau=1", "[traj][motion_laws]")
{
    auto pt = septic_law(1.0);
    CHECK_THAT(pt.q, WithinAbs(1.0, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(0.0, 1e-12));
}

TEST_CASE("septic_law at tau=0.5 (symmetric midpoint)", "[traj][motion_laws]")
{
    auto pt = septic_law(0.5);
    CHECK_THAT(pt.q, WithinAbs(0.5, 1e-12));
}

// --- Harmonic law ---

TEST_CASE("harmonic_law at tau=0", "[traj][motion_laws]")
{
    auto const pi2 = std::numbers::pi * std::numbers::pi;
    auto pt = harmonic_law(0.0);
    CHECK_THAT(pt.q, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(pi2 / 2.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(0.0, 1e-12));
}

TEST_CASE("harmonic_law at tau=1", "[traj][motion_laws]")
{
    auto const pi2 = std::numbers::pi * std::numbers::pi;
    auto pt = harmonic_law(1.0);
    CHECK_THAT(pt.q, WithinAbs(1.0, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(-pi2 / 2.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(0.0, 1e-12));
}

TEST_CASE("harmonic_law at tau=0.5 (symmetric midpoint)", "[traj][motion_laws]")
{
    auto pt = harmonic_law(0.5);
    CHECK_THAT(pt.q, WithinAbs(0.5, 1e-12));
}

// --- Cycloidal law ---

TEST_CASE("cycloidal_law at tau=0", "[traj][motion_laws]")
{
    auto const pi2 = std::numbers::pi * std::numbers::pi;
    auto pt = cycloidal_law(0.0);
    CHECK_THAT(pt.q, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(4.0 * pi2, 1e-12));
}

TEST_CASE("cycloidal_law at tau=1", "[traj][motion_laws]")
{
    auto const pi2 = std::numbers::pi * std::numbers::pi;
    auto pt = cycloidal_law(1.0);
    CHECK_THAT(pt.q, WithinAbs(1.0, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(4.0 * pi2, 1e-12));
}

TEST_CASE("cycloidal_law at tau=0.5 (symmetric midpoint)", "[traj][motion_laws]")
{
    auto pt = cycloidal_law(0.5);
    CHECK_THAT(pt.q, WithinAbs(0.5, 1e-12));
}

// --- Monotonicity for all laws ---

TEST_CASE("all motion laws are monotonically increasing", "[traj][motion_laws]")
{
    auto check_monotonic = [](auto law, const char* name) {
        auto q25 = law(0.25).q;
        auto q50 = law(0.5).q;
        auto q75 = law(0.75).q;
        INFO("Law: " << name);
        CHECK(q25 < q50);
        CHECK(q50 < q75);
    };

    check_monotonic(cubic_law<double>, "cubic");
    check_monotonic(quintic_law<double>, "quintic");
    check_monotonic(septic_law<double>, "septic");
    check_monotonic(harmonic_law<double>, "harmonic");
    check_monotonic(cycloidal_law<double>, "cycloidal");
}
