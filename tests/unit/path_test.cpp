#include "ctrlpp/traj/cubic_path.h"
#include "ctrlpp/traj/cycloidal_path.h"
#include "ctrlpp/traj/harmonic_path.h"
#include "ctrlpp/traj/quintic_path.h"
#include "ctrlpp/traj/septic_path.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <numbers>

using namespace ctrlpp;
using Catch::Matchers::WithinAbs;

// --- Cubic path ---

TEST_CASE("cubic_path at tau=0", "[traj][path]")
{
    auto pt = cubic_path(0.0);
    CHECK_THAT(pt.q, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(6.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(-12.0, 1e-12));
}

TEST_CASE("cubic_path at tau=1", "[traj][path]")
{
    auto pt = cubic_path(1.0);
    CHECK_THAT(pt.q, WithinAbs(1.0, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(-6.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(-12.0, 1e-12));
}

TEST_CASE("cubic_path at tau=0.5", "[traj][path]")
{
    auto pt = cubic_path(0.5);
    CHECK_THAT(pt.q, WithinAbs(0.5, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(1.5, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(-12.0, 1e-12));
}

// --- Quintic path ---

TEST_CASE("quintic_path at tau=0", "[traj][path]")
{
    auto pt = quintic_path(0.0);
    CHECK_THAT(pt.q, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(60.0, 1e-12));
}

TEST_CASE("quintic_path at tau=1", "[traj][path]")
{
    auto pt = quintic_path(1.0);
    CHECK_THAT(pt.q, WithinAbs(1.0, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(60.0, 1e-12));
}

TEST_CASE("quintic_path at tau=0.5 (symmetric midpoint)", "[traj][path]")
{
    auto pt = quintic_path(0.5);
    CHECK_THAT(pt.q, WithinAbs(0.5, 1e-12));
}

// --- Septic path ---

TEST_CASE("septic_path at tau=0", "[traj][path]")
{
    auto pt = septic_path(0.0);
    CHECK_THAT(pt.q, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(0.0, 1e-12));
}

TEST_CASE("septic_path at tau=1", "[traj][path]")
{
    auto pt = septic_path(1.0);
    CHECK_THAT(pt.q, WithinAbs(1.0, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(0.0, 1e-12));
}

TEST_CASE("septic_path at tau=0.5 (symmetric midpoint)", "[traj][path]")
{
    auto pt = septic_path(0.5);
    CHECK_THAT(pt.q, WithinAbs(0.5, 1e-12));
}

// --- Harmonic path ---

TEST_CASE("harmonic_path at tau=0", "[traj][path]")
{
    auto const pi2 = std::numbers::pi * std::numbers::pi;
    auto pt = harmonic_path(0.0);
    CHECK_THAT(pt.q, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(pi2 / 2.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(0.0, 1e-12));
}

TEST_CASE("harmonic_path at tau=1", "[traj][path]")
{
    auto const pi2 = std::numbers::pi * std::numbers::pi;
    auto pt = harmonic_path(1.0);
    CHECK_THAT(pt.q, WithinAbs(1.0, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(-pi2 / 2.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(0.0, 1e-12));
}

TEST_CASE("harmonic_path at tau=0.5 (symmetric midpoint)", "[traj][path]")
{
    auto pt = harmonic_path(0.5);
    CHECK_THAT(pt.q, WithinAbs(0.5, 1e-12));
}

// --- Cycloidal path ---

TEST_CASE("cycloidal_path at tau=0", "[traj][path]")
{
    auto const pi2 = std::numbers::pi * std::numbers::pi;
    auto pt = cycloidal_path(0.0);
    CHECK_THAT(pt.q, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(4.0 * pi2, 1e-12));
}

TEST_CASE("cycloidal_path at tau=1", "[traj][path]")
{
    auto const pi2 = std::numbers::pi * std::numbers::pi;
    auto pt = cycloidal_path(1.0);
    CHECK_THAT(pt.q, WithinAbs(1.0, 1e-12));
    CHECK_THAT(pt.dq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.ddq, WithinAbs(0.0, 1e-12));
    CHECK_THAT(pt.dddq, WithinAbs(4.0 * pi2, 1e-12));
}

TEST_CASE("cycloidal_path at tau=0.5 (symmetric midpoint)", "[traj][path]")
{
    auto pt = cycloidal_path(0.5);
    CHECK_THAT(pt.q, WithinAbs(0.5, 1e-12));
}

// --- Monotonicity for all paths ---

TEST_CASE("all paths are monotonically increasing", "[traj][path]")
{
    auto check_monotonic = [](auto path, const char* name) {
        auto q25 = path(0.25).q;
        auto q50 = path(0.5).q;
        auto q75 = path(0.75).q;
        INFO("Path: " << name);
        CHECK(q25 < q50);
        CHECK(q50 < q75);
    };

    check_monotonic(cubic_path<double>, "cubic");
    check_monotonic(quintic_path<double>, "quintic");
    check_monotonic(septic_path<double>, "septic");
    check_monotonic(harmonic_path<double>, "harmonic");
    check_monotonic(cycloidal_path<double>, "cycloidal");
}

// --- Peak derivative tests ---

TEST_CASE("cubic_path_peak_derivatives returns expected values", "[traj][path][peaks]")
{
    auto peaks = cubic_path_peak_derivatives<double>();
    CHECK_THAT(peaks[0], WithinAbs(1.5, 1e-12));
    CHECK_THAT(peaks[1], WithinAbs(6.0, 1e-12));
    CHECK_THAT(peaks[2], WithinAbs(12.0, 1e-12));
}

TEST_CASE("quintic_path_peak_derivatives returns expected values", "[traj][path][peaks]")
{
    auto peaks = quintic_path_peak_derivatives<double>();
    CHECK_THAT(peaks[0], WithinAbs(15.0 / 8.0, 1e-12));
    CHECK_THAT(peaks[1], WithinAbs(10.0 * std::sqrt(3.0) / 3.0, 1e-10));
    CHECK_THAT(peaks[2], WithinAbs(60.0, 1e-12));
}

TEST_CASE("septic_path_peak_derivatives returns expected values", "[traj][path][peaks]")
{
    auto peaks = septic_path_peak_derivatives<double>();
    CHECK_THAT(peaks[0], WithinAbs(35.0 / 16.0, 1e-12));
    CHECK_THAT(peaks[1], WithinAbs(7.5132, 1e-3));
    CHECK_THAT(peaks[2], WithinAbs(52.5, 1e-12));
}

TEST_CASE("harmonic_path_peak_derivatives returns expected values", "[traj][path][peaks]")
{
    auto const pi = std::numbers::pi;
    auto peaks = harmonic_path_peak_derivatives<double>();
    CHECK_THAT(peaks[0], WithinAbs(pi / 2.0, 1e-12));
    CHECK_THAT(peaks[1], WithinAbs(pi * pi / 2.0, 1e-12));
    CHECK_THAT(peaks[2], WithinAbs(pi * pi * pi / 2.0, 1e-12));
}

TEST_CASE("cycloidal_path_peak_derivatives returns expected values", "[traj][path][peaks]")
{
    auto const pi = std::numbers::pi;
    auto peaks = cycloidal_path_peak_derivatives<double>();
    CHECK_THAT(peaks[0], WithinAbs(2.0, 1e-12));
    CHECK_THAT(peaks[1], WithinAbs(2.0 * pi, 1e-12));
    CHECK_THAT(peaks[2], WithinAbs(4.0 * pi * pi, 1e-12));
}
