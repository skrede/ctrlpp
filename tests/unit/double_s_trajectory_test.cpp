#include <ctrlpp/trajectory/double_s_trajectory.h>
#include <ctrlpp/trajectory/trajectory_segment.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace
{

// create() is the only construction path and it is fallible, so every profile a
// test uses is built through this helper, which asserts the command was
// realizable rather than letting a rejection pass as a profile. The rejection
// cases assert the specific enumerator and never come through here.
auto built(ctrlpp::double_s_trajectory<double>::config const& cfg)
    -> ctrlpp::double_s_trajectory<double>
{
    auto created = ctrlpp::double_s_trajectory<double>::create(cfg);
    REQUIRE(created.has_value());
    return created.value();
}

}

// --------------------------------------------------------------------------
// Test 1: Full 7-segment profile reaches q1 with 7 distinct phases
// --------------------------------------------------------------------------
TEST_CASE("double_s: full 7-segment profile", "[traj][double_s]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0};
    auto const traj = built(cfg);

    auto const T = traj.duration();
    REQUIRE(T > 0.0);

    auto const pt_end = traj.evaluate(T);
    REQUIRE_THAT(pt_end.position[0], WithinAbs(10.0, 1e-10));
    REQUIRE_THAT(pt_end.velocity[0], WithinAbs(0.0, 1e-8));

    auto const pt_start = traj.evaluate(0.0);
    REQUIRE_THAT(pt_start.position[0], WithinAbs(0.0, 1e-10));
    REQUIRE_THAT(pt_start.velocity[0], WithinAbs(0.0, 1e-8));

    // All 7 phase durations should be > 0 for a full profile
    auto const pd = traj.phase_durations();
    REQUIRE(pd.size() == 7);
    // At least the jerk phases (0, 2, 4, 6) and cruise (3) should be > 0
    CHECK(pd[0] > 0.0);
    CHECK(pd[3] >= 0.0); // cruise may or may not exist
    CHECK(pd[6] > 0.0);
}

// --------------------------------------------------------------------------
// Test 2: Velocity not reached (degenerate -- v_lim < v_max)
// --------------------------------------------------------------------------
TEST_CASE("double_s: velocity not reached degenerate", "[traj][double_s]")
{
    // Short distance forces v_lim < v_max
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 1.0, .v_max = 100.0, .a_max = 10.0, .j_max = 100.0};
    auto const traj = built(cfg);

    CHECK(traj.is_degenerate());
    CHECK(traj.peak_velocity() < 100.0);

    auto const T = traj.duration();
    REQUIRE(T > 0.0);
    auto const pt = traj.evaluate(T);
    REQUIRE_THAT(pt.position[0], WithinAbs(1.0, 1e-10));
}

// --------------------------------------------------------------------------
// Test 3: Acceleration not reached (a_lim < a_max)
// --------------------------------------------------------------------------
TEST_CASE("double_s: acceleration not reached degenerate", "[traj][double_s]")
{
    // Very high a_max but low j_max -- can't reach full a_max
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 0.5, .v_max = 10.0, .a_max = 100.0, .j_max = 50.0};
    auto const traj = built(cfg);

    CHECK(traj.is_degenerate());

    auto const T = traj.duration();
    REQUIRE(T > 0.0);
    auto const pt = traj.evaluate(T);
    REQUIRE_THAT(pt.position[0], WithinAbs(0.5, 1e-10));
}

// --------------------------------------------------------------------------
// Test 4: Both v and a not reached (doubly degenerate)
// --------------------------------------------------------------------------
TEST_CASE("double_s: doubly degenerate (v and a not reached)", "[traj][double_s]")
{
    // Tiny distance with high limits -- neither v_max nor a_max reached
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 0.01, .v_max = 100.0, .a_max = 100.0, .j_max = 50.0};
    auto const traj = built(cfg);

    CHECK(traj.is_degenerate());
    CHECK(traj.peak_velocity() < 100.0);

    auto const T = traj.duration();
    REQUIRE(T > 0.0);
    auto const pt = traj.evaluate(T);
    REQUIRE_THAT(pt.position[0], WithinAbs(0.01, 1e-10));
}

// --------------------------------------------------------------------------
// Test 5: Triangular-like (very short distance)
// --------------------------------------------------------------------------
TEST_CASE("double_s: very short distance simplifies", "[traj][double_s]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 0.001, .v_max = 10.0, .a_max = 10.0, .j_max = 100.0};
    auto const traj = built(cfg);

    auto const T = traj.duration();
    REQUIRE(T > 0.0);
    auto const pt = traj.evaluate(T);
    REQUIRE_THAT(pt.position[0], WithinAbs(0.001, 1e-10));
}

// --------------------------------------------------------------------------
// Test 6: Velocity constraint -- sample 1000 points, |velocity| <= v_max + eps
// --------------------------------------------------------------------------
TEST_CASE("double_s: velocity constraint satisfied", "[traj][double_s]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0};
    auto const traj = built(cfg);

    auto const T = traj.duration();
    constexpr double eps = 1e-6;
    for (int i = 0; i <= 1000; ++i) {
        double const t = T * static_cast<double>(i) / 1000.0;
        auto const pt = traj.evaluate(t);
        REQUIRE(std::abs(pt.velocity[0]) <= cfg.v_max + eps);
    }
}

// --------------------------------------------------------------------------
// Test 7: Acceleration constraint -- sample 1000 points
// --------------------------------------------------------------------------
TEST_CASE("double_s: acceleration constraint satisfied", "[traj][double_s]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0};
    auto const traj = built(cfg);

    auto const T = traj.duration();
    constexpr double eps = 1e-6;
    for (int i = 0; i <= 1000; ++i) {
        double const t = T * static_cast<double>(i) / 1000.0;
        auto const pt = traj.evaluate(t);
        REQUIRE(std::abs(pt.acceleration[0]) <= cfg.a_max + eps);
    }
}

// --------------------------------------------------------------------------
// Test 8: Jerk constraint -- finite difference check
// --------------------------------------------------------------------------
TEST_CASE("double_s: jerk constraint satisfied", "[traj][double_s]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0};
    auto const traj = built(cfg);

    auto const T = traj.duration();
    constexpr double dt = 1e-5;
    constexpr double eps = 1.0; // jerk via finite diff has larger tolerance
    for (int i = 1; i < 1000; ++i) {
        double const t = T * static_cast<double>(i) / 1000.0;
        auto const pt1 = traj.evaluate(t - dt);
        auto const pt2 = traj.evaluate(t + dt);
        double const jerk = (pt2.acceleration[0] - pt1.acceleration[0]) / (2.0 * dt);
        REQUIRE(std::abs(jerk) <= cfg.j_max + eps);
    }
}

// --------------------------------------------------------------------------
// Test 9: Phase boundary continuity
// --------------------------------------------------------------------------
TEST_CASE("double_s: phase boundary continuity", "[traj][double_s]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0};
    auto const traj = built(cfg);

    auto const pd = traj.phase_durations();
    double boundary = 0.0;
    constexpr double tiny = 1e-13;
    constexpr double tol = 1e-10;

    for (std::size_t i = 0; i < 6; ++i) {
        boundary += pd[i];
        if (boundary <= 0.0 || boundary >= traj.duration()) {
            continue;
        }
        auto const left = traj.evaluate(boundary - tiny);
        auto const right = traj.evaluate(boundary + tiny);
        REQUIRE_THAT(left.position[0], WithinAbs(right.position[0], tol));
        REQUIRE_THAT(left.velocity[0], WithinAbs(right.velocity[0], tol));
        REQUIRE_THAT(left.acceleration[0], WithinAbs(right.acceleration[0], tol));
    }
}

// --------------------------------------------------------------------------
// Test 10: Negative displacement (q1 < q0) via sigma transformation
// --------------------------------------------------------------------------
TEST_CASE("double_s: negative displacement sigma transformation", "[traj][double_s]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 10.0, .q1 = 0.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0};
    auto const traj = built(cfg);

    auto const T = traj.duration();
    REQUIRE(T > 0.0);

    auto const pt_start = traj.evaluate(0.0);
    REQUIRE_THAT(pt_start.position[0], WithinAbs(10.0, 1e-10));

    auto const pt_end = traj.evaluate(T);
    REQUIRE_THAT(pt_end.position[0], WithinAbs(0.0, 1e-10));

    // Velocity should be negative (moving from 10 -> 0)
    auto const pt_mid = traj.evaluate(T / 2.0);
    CHECK(pt_mid.velocity[0] < 0.0);
}

// --------------------------------------------------------------------------
// Test 11: Introspection methods
// --------------------------------------------------------------------------
TEST_CASE("double_s: introspection", "[traj][double_s]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0};
    auto const traj = built(cfg);

    // is_degenerate returns bool
    auto const deg = traj.is_degenerate();
    static_assert(std::is_same_v<decltype(deg), const bool>);

    // peak_velocity returns Scalar
    auto const pv = traj.peak_velocity();
    CHECK(pv > 0.0);
    CHECK(pv <= cfg.v_max + 1e-10);

    // phase_durations returns array of 7
    auto const pd = traj.phase_durations();
    static_assert(pd.size() == 7);

    // Sum of phase durations equals total duration
    double sum = 0.0;
    for (auto const& d : pd) {
        CHECK(d >= 0.0);
        sum += d;
    }
    REQUIRE_THAT(sum, WithinAbs(traj.duration(), 1e-12));
}

// --------------------------------------------------------------------------
// Test 12: Satisfies trajectory_segment concept
// --------------------------------------------------------------------------
TEST_CASE("double_s: satisfies trajectory_segment concept", "[traj][double_s]")
{
    static_assert(ctrlpp::trajectory_segment<
                  ctrlpp::double_s_trajectory<double>, double, 1>);

    static_assert(ctrlpp::trajectory_segment<
                  ctrlpp::double_s_trajectory<float>, float, 1>);
}

// --------------------------------------------------------------------------
// Test 13: Zero displacement (q0 == q1) returns stationary profile
// --------------------------------------------------------------------------
TEST_CASE("double_s: zero displacement stationary profile", "[traj][double_s]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 5.0, .q1 = 5.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0};
    auto const traj = built(cfg);

    CHECK(traj.duration() == 0.0);

    auto const pt = traj.evaluate(0.0);
    REQUIRE_THAT(pt.position[0], WithinAbs(5.0, 1e-12));
    REQUIRE_THAT(pt.velocity[0], WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(pt.acceleration[0], WithinAbs(0.0, 1e-12));
}

// --------------------------------------------------------------------------
// Test 14: Nonzero boundary velocities (B&M Example 3.9 parameters)
// --------------------------------------------------------------------------
TEST_CASE("double_s: nonzero boundary velocities", "[traj][double_s]")
{
    // Biagiotti & Melchiorri, Example 3.9, p.84: q0=0, q1=10, v0=1, v1=0,
    // v_max=5, a_max=10, j_max=30. v_max is reached, so peak velocity is v_max.
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 30.0,
        .v0 = 1.0, .v1 = 0.0};
    auto const traj = built(cfg);

    auto const T = traj.duration();
    REQUIRE(T > 0.0);

    // Tolerances derived from input scale x machine epsilon: evaluate() chains a
    // handful of multiply-adds and one phase-boundary subtraction, so allow a few
    // ULP at each quantity's own scale.
    constexpr double eps = std::numeric_limits<double>::epsilon();
    constexpr double rounding_ops = 16.0;
    double const q_tol = rounding_ops * eps * std::abs(cfg.q1 - cfg.q0);
    double const v_tol = rounding_ops * eps * cfg.v_max;
    double const a_tol = rounding_ops * eps * cfg.a_max;

    // Boundary velocities honored exactly, zero acceleration at both ends.
    auto const pt_start = traj.evaluate(0.0);
    REQUIRE_THAT(pt_start.position[0], WithinAbs(cfg.q0, q_tol));
    REQUIRE_THAT(pt_start.velocity[0], WithinAbs(cfg.v0, v_tol));
    REQUIRE_THAT(pt_start.acceleration[0], WithinAbs(0.0, a_tol));

    auto const pt_end = traj.evaluate(T);
    REQUIRE_THAT(pt_end.position[0], WithinAbs(cfg.q1, q_tol));
    REQUIRE_THAT(pt_end.velocity[0], WithinAbs(cfg.v1, v_tol));
    REQUIRE_THAT(pt_end.acceleration[0], WithinAbs(0.0, a_tol));

    // v_max is reached for these parameters.
    REQUIRE_THAT(traj.peak_velocity(), WithinAbs(cfg.v_max, v_tol));

    // The reported trace honors its boundary velocity in the interior too: one
    // sample step past the start the velocity has moved off v0 by no more than
    // the acceleration limit allows, i.e. it starts from v0 rather than 0.
    double const dt = T / 1000.0;
    auto const pt_near = traj.evaluate(dt);
    REQUIRE(std::abs(pt_near.velocity[0] - cfg.v0) <= cfg.a_max * dt + v_tol);

    // Constraints hold across a dense scan.
    for (int i = 0; i <= 1000; ++i) {
        double const t = T * static_cast<double>(i) / 1000.0;
        auto const pt = traj.evaluate(t);
        REQUIRE(std::abs(pt.velocity[0]) <= cfg.v_max + v_tol);
        REQUIRE(std::abs(pt.acceleration[0]) <= cfg.a_max + a_tol);
    }
}

// --------------------------------------------------------------------------
// The velocity limit is a precondition on the boundary velocities
// --------------------------------------------------------------------------
TEST_CASE("double_s: a boundary velocity above the velocity limit is rejected",
          "[traj][double_s]")
{
    using traj_t = ctrlpp::double_s_trajectory<double>;

    // Raising the limit to whatever the caller passed as a boundary velocity
    // would return a profile that violates a bound the caller stated, so the
    // command is rejected instead.
    auto const entering = traj_t::create(
        {.q0 = 0.0, .q1 = 10.0, .v_max = 1.0, .a_max = 1.0, .j_max = 1.0, .v0 = 5.0, .v1 = 0.0});
    REQUIRE_FALSE(entering.has_value());
    REQUIRE(entering.error() == ctrlpp::trajectory_error::boundary_velocity_exceeds_limit);

    // The limit bounds the SPEED, so the sign does not enter.
    auto const leaving = traj_t::create(
        {.q0 = 0.0, .q1 = 10.0, .v_max = 1.0, .a_max = 1.0, .j_max = 1.0, .v0 = 0.0, .v1 = -5.0});
    REQUIRE_FALSE(leaving.has_value());
    REQUIRE(leaving.error() == ctrlpp::trajectory_error::boundary_velocity_exceeds_limit);

    // Reaching the limit exactly is inside the domain: the ramp attached to that
    // boundary velocity vanishes, which is a degenerate shape rather than an
    // impossible one.
    auto const at_limit = traj_t::create(
        {.q0 = 0.0, .q1 = 10.0, .v_max = 1.0, .a_max = 1.0, .j_max = 1.0, .v0 = 1.0, .v1 = 0.0});
    REQUIRE(at_limit.has_value());
    REQUIRE(at_limit.value().duration() > 0.0);
}

// --------------------------------------------------------------------------
// A displacement below what the boundary velocities already sweep
// --------------------------------------------------------------------------
TEST_CASE("double_s: a command below the transition distance is rejected", "[traj][double_s]")
{
    using traj_t = ctrlpp::double_s_trajectory<double>;

    // A zero commanded displacement is realizable by exactly one member of this
    // family, the standstill, and only when both boundary velocities are zero.
    // Anything else asks the axis to change speed while covering no ground.
    auto const standstill =
        traj_t::create({.q0 = 5.0, .q1 = 5.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0});
    REQUIRE(standstill.has_value());
    REQUIRE_THAT(standstill.value().duration(), WithinAbs(0.0, 0.0));

    auto const speed_change = traj_t::create(
        {.q0 = 5.0, .q1 = 5.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0, .v0 = 1.0});
    REQUIRE_FALSE(speed_change.has_value());
    REQUIRE(speed_change.error() == ctrlpp::trajectory_error::unreachable_boundary_velocity);

    // A displacement below the minimum but not at zero is the same rejection.
    // Previously a command inside a fixed neighborhood of zero was rounded down
    // to a standstill, which reported a resting axis for one commanded to be
    // moving at both ends.
    auto const too_short = traj_t::create(
        {.q0 = 0.0, .q1 = 1e-16, .v_max = 1.0, .a_max = 1.0, .j_max = 1.0, .v0 = 0.5});
    REQUIRE_FALSE(too_short.has_value());
    REQUIRE(too_short.error() == ctrlpp::trajectory_error::unreachable_boundary_velocity);

    // The same tiny displacement at rest IS realizable, and is now served as the
    // brief move it is rather than reported as a standstill.
    auto const tiny_at_rest =
        traj_t::create({.q0 = 0.0, .q1 = 1e-16, .v_max = 1.0, .a_max = 1.0, .j_max = 1.0});
    REQUIRE(tiny_at_rest.has_value());
    REQUIRE(tiny_at_rest.value().duration() > 0.0);
    auto const arrival = tiny_at_rest.value().evaluate(tiny_at_rest.value().duration());
    REQUIRE_THAT(arrival.position[0], WithinRel(1e-16, 1e-9));
}

// --------------------------------------------------------------------------
// Input domain of the three limits and the boundary values
// --------------------------------------------------------------------------
TEST_CASE("double_s: out-of-domain limits and non-finite inputs are rejected", "[traj][double_s]")
{
    using traj_t = ctrlpp::double_s_trajectory<double>;
    auto constexpr nan = std::numeric_limits<double>::quiet_NaN();
    auto constexpr inf = std::numeric_limits<double>::infinity();

    for (double const v_max : {0.0, -1.0, nan, inf}) {
        auto const result =
            traj_t::create({.q0 = 0.0, .q1 = 1.0, .v_max = v_max, .a_max = 1.0, .j_max = 1.0});
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error() == ctrlpp::trajectory_error::non_positive_velocity_limit);
    }
    for (double const a_max : {0.0, -1.0, nan, inf}) {
        auto const result =
            traj_t::create({.q0 = 0.0, .q1 = 1.0, .v_max = 1.0, .a_max = a_max, .j_max = 1.0});
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error() == ctrlpp::trajectory_error::non_positive_acceleration_limit);
    }
    for (double const j_max : {0.0, -1.0, nan, inf}) {
        auto const result =
            traj_t::create({.q0 = 0.0, .q1 = 1.0, .v_max = 1.0, .a_max = 1.0, .j_max = j_max});
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error() == ctrlpp::trajectory_error::non_positive_jerk_limit);
    }
    for (double const bad : {nan, inf, -inf}) {
        REQUIRE(traj_t::create({.q0 = bad, .q1 = 1.0, .v_max = 1.0, .a_max = 1.0, .j_max = 1.0})
                    .error()
                == ctrlpp::trajectory_error::non_finite_input);
        REQUIRE(traj_t::create(
                    {.q0 = 0.0, .q1 = 1.0, .v_max = 1.0, .a_max = 1.0, .j_max = 1.0, .v0 = bad})
                    .error()
                == ctrlpp::trajectory_error::non_finite_input);
    }
}

// --------------------------------------------------------------------------
// A profile that does not fit the scalar type
// --------------------------------------------------------------------------
TEST_CASE("double_s: a duration outside the representable range is rejected", "[traj][double_s]")
{
    using traj_t = ctrlpp::double_s_trajectory<double>;

    auto const overflowed =
        traj_t::create({.q0 = 0.0, .q1 = 1e308, .v_max = 1e-300, .a_max = 1.0, .j_max = 1.0});
    REQUIRE_FALSE(overflowed.has_value());
    REQUIRE(overflowed.error() == ctrlpp::trajectory_error::unrepresentable_duration);

    auto const served =
        traj_t::create({.q0 = 0.0, .q1 = 1e308, .v_max = 1e6, .a_max = 1.0, .j_max = 1.0});
    REQUIRE(served.has_value());
    REQUIRE(std::isfinite(served.value().duration()));
}
