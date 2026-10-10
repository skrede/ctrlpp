#include <ctrlpp/trajectory/double_s_trajectory.h>
#include <ctrlpp/trajectory/trajectory_segment.h>

#include <ctrlpp/trajectory/detail/jerk_limited_ramp.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bit>
#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <cstdint>
#include <algorithm>
#include <type_traits>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinULP;

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

template <typename Scalar>
using scalar_bits = std::conditional_t<sizeof(Scalar) == sizeof(std::uint64_t), std::uint64_t, std::uint32_t>;

template <typename Scalar>
auto bits_of(Scalar value) -> scalar_bits<Scalar>
{
    return std::bit_cast<scalar_bits<Scalar>>(value);
}

// The out-parameters start as NaN so a component the overload never writes
// cannot pass by coinciding with the by-value result.
template <typename Scalar>
void require_bit_identical_at(ctrlpp::double_s_trajectory<Scalar> const& traj, Scalar t)
{
    auto const by_value = traj.evaluate(t);
    auto q = std::numeric_limits<Scalar>::quiet_NaN();
    auto dq = q;
    auto ddq = q;
    traj.evaluate(t, q, dq, ddq);

    CAPTURE(t);
    REQUIRE(bits_of(q) == bits_of(by_value.position[0]));
    REQUIRE(bits_of(dq) == bits_of(by_value.velocity[0]));
    REQUIRE(bits_of(ddq) == bits_of(by_value.acceleration[0]));
}

// The sweep runs from an eighth of the duration before the start to an eighth
// past the end, so both clamped tails are sampled along with every segment.
template <typename Scalar>
void require_bit_identical_over_horizon(ctrlpp::double_s_trajectory<Scalar> const& traj)
{
    constexpr std::int32_t samples = 1025;
    auto const T = traj.duration();
    for (std::int32_t i = 0; i < samples; ++i) {
        auto const fraction = static_cast<Scalar>(i) / static_cast<Scalar>(samples - 1);
        require_bit_identical_at(traj, -T / Scalar{8} + T * Scalar{1.25} * fraction);
    }
    for (Scalar const t : {Scalar{0}, T, -T, Scalar{2} * T, Scalar{-1}, Scalar{1}}) {
        require_bit_identical_at(traj, t);
    }
}

template <typename Scalar>
auto built_profile(typename ctrlpp::double_s_trajectory<Scalar>::config const& cfg)
    -> ctrlpp::double_s_trajectory<Scalar>
{
    auto created = ctrlpp::double_s_trajectory<Scalar>::create(cfg);
    REQUIRE(created.has_value());
    return created.value();
}

auto sampled_peak_speed(ctrlpp::double_s_trajectory<double> const& traj) -> double
{
    auto const T = traj.duration();
    auto peak = std::abs(traj.evaluate(T).velocity[0]);
    for (std::int64_t i = 0; static_cast<double>(i) * 1e-3 < T; ++i) {
        peak = std::max(peak, std::abs(traj.evaluate(static_cast<double>(i) * 1e-3).velocity[0]));
    }
    return peak;
}

struct limit_set
{
    double v_max, a_max, j_max;
};

constexpr std::array<limit_set, 4> near_limit_sets{{{1.0, 1.0, 1.0}, {5.0, 10.0, 50.0}, {0.7, 1.3, 2.9}, {1.0, 5.0, 50.0}}};

auto stepped(double value, std::int32_t steps) -> double
{
    auto const toward = (steps < 0) ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity();
    for (std::int32_t i = 0; i < std::abs(steps); ++i) {
        value = std::nextafter(value, toward);
    }
    return value;
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
    REQUIRE(at_limit.value().peak_velocity() <= 1.0);
    // Two roundings: the subtraction of the swept distance and the division by v_max.
    REQUIRE_THAT(at_limit.value().phase_durations()[3], WithinULP(9.0, 2));
}

TEST_CASE("Double-S: a boundary velocity at the limit builds a profile at the limit",
          "[trajectory][double_s]")
{
    struct frame
    {
        double q1, v0, v1;
    };
    for (auto const f : std::array<frame, 4>{{{10.0, 1.0, 0.0}, {10.0, 0.0, 1.0}, {-10.0, -1.0, 0.0}, {-10.0, 0.0, -1.0}}}) {
        CAPTURE(f.q1, f.v0, f.v1);
        auto const traj =
            built({.q0 = 0.0, .q1 = f.q1, .v_max = 1.0, .a_max = 1.0, .j_max = 1.0, .v0 = f.v0, .v1 = f.v1});
        REQUIRE(traj.peak_velocity() <= 1.0);
        // Four roundings: the subtraction of the swept distance, the division by
        // v_max and the two sums that total the segments.
        REQUIRE_THAT(traj.duration(), WithinULP(11.0, 4));
        // Two roundings: the subtraction of the swept distance and the division by v_max.
        REQUIRE_THAT(traj.phase_durations()[3], WithinULP(9.0, 2));
        REQUIRE(sampled_peak_speed(traj) <= 1.0);
    }
}

TEST_CASE("Double-S: a rest-to-rest peak that rounds above the limit is held at the limit", "[trajectory][double_s]")
{
    // Two representable steps below the cruise boundary; v_max * j_max == a_max^2
    // sends it to the doubly degenerate branch.
    auto const traj = built({.q0 = 0.0, .q1 = 1.9999999999999996, .v_max = 1.0, .a_max = 1.0, .j_max = 1.0});
    REQUIRE(traj.peak_velocity() <= 1.0);
}

TEST_CASE("Double-S: rest-to-rest commands around the cruise boundary build at or below the limit",
          "[trajectory][double_s]")
{
    for (auto const l : near_limit_sets) {
        // @cite biagiotti2009 -- Sec. 3.4.3, p.89
        auto const T_a = (l.v_max * l.j_max >= l.a_max * l.a_max) ? l.a_max / l.j_max + l.v_max / l.a_max
                                                                   : 2.0 * std::sqrt(l.v_max / l.j_max);
        auto h = stepped(l.v_max * T_a, -400);
        for (std::int32_t i = 0; i <= 800; ++i, h = stepped(h, 1)) {
            CAPTURE(l.v_max, l.a_max, l.j_max, h);
            auto const created = ctrlpp::double_s_trajectory<double>::create(
                {.q0 = 0.0, .q1 = h, .v_max = l.v_max, .a_max = l.a_max, .j_max = l.j_max});
            REQUIRE(created.has_value());
            REQUIRE(created.value().peak_velocity() <= l.v_max);
        }
    }
}

TEST_CASE("Double-S: entering at the limit near the fastest transition builds at or below the limit or is unreachable",
          "[trajectory][double_s]")
{
    for (auto const l : near_limit_sets) {
        for (double const fraction : {0.0, 0.25, 0.5, 0.75, 1.0}) {
            auto const v1 = fraction * l.v_max;
            auto const dv = l.v_max - v1;
            // @cite biagiotti2009 -- Sec. 3.4.1, eq. (3.19)-(3.20), p.80
            auto const T = (dv * l.j_max < l.a_max * l.a_max) ? 2.0 * std::sqrt(dv / l.j_max) : dv / l.a_max + l.a_max / l.j_max;
            auto h = stepped((l.v_max + v1) * T / 2.0, -20);
            for (std::int32_t i = 0; i <= 420; ++i, h = stepped(h, 1)) {
                CAPTURE(l.v_max, l.a_max, l.j_max, v1, h);
                auto const created = ctrlpp::double_s_trajectory<double>::create(
                    {.q0 = 0.0, .q1 = h, .v_max = l.v_max, .a_max = l.a_max, .j_max = l.j_max, .v0 = l.v_max, .v1 = v1});
                if (created.has_value()) {
                    REQUIRE(created.value().peak_velocity() <= l.v_max);
                } else if (h != 0.0 && h / l.v_max == 0.0) {
                    // With both boundary velocities at the limit the fastest transition is zero, so the scan
                    // reaches subnormal displacements whose traversal time at the limit underflows.
                    REQUIRE(created.error() == ctrlpp::trajectory_error::unrepresentable_duration);
                } else {
                    REQUIRE(created.error() == ctrlpp::trajectory_error::unreachable_boundary_velocity);
                }
            }
        }
    }
}

TEMPLATE_TEST_CASE("Double-S: a counted boundary velocity tolerance admits and clamps, and refuses beyond it",
                   "[trajectory][double_s]", float, double)
{
    using Scalar = TestType;
    using traj_t = ctrlpp::double_s_trajectory<Scalar>;
    auto const create = [](Scalar q1, Scalar v0, Scalar v1, std::uint32_t rounding_ops) {
        return traj_t::create({.q0 = Scalar{0},
                               .q1 = q1,
                               .v_max = Scalar{1},
                               .a_max = Scalar{1},
                               .j_max = Scalar{1},
                               .v0 = v0,
                               .v1 = v1,
                               .boundary_velocity_rounding_ops = rounding_ops});
    };
    auto const one_ulp_above = std::nextafter(Scalar{1}, Scalar{2});
    auto const two_ulp_above = std::nextafter(one_ulp_above, Scalar{2});

    auto const exact = create(Scalar{10}, one_ulp_above, Scalar{0}, 0);
    REQUIRE_FALSE(exact.has_value());
    REQUIRE(exact.error() == ctrlpp::trajectory_error::boundary_velocity_exceeds_limit);

    auto const entering = create(Scalar{10}, one_ulp_above, Scalar{0}, 1);
    REQUIRE(entering.has_value());
    REQUIRE(entering.value().evaluate(Scalar{0}).velocity[0] == Scalar{1});
    REQUIRE(entering.value().peak_velocity() <= Scalar{1});

    auto const beyond = create(Scalar{10}, two_ulp_above, Scalar{0}, 1);
    REQUIRE_FALSE(beyond.has_value());
    REQUIRE(beyond.error() == ctrlpp::trajectory_error::boundary_velocity_exceeds_limit);
    REQUIRE(create(Scalar{10}, two_ulp_above, Scalar{0}, 2).has_value());

    auto const negative = create(Scalar{-10}, -one_ulp_above, Scalar{0}, 1);
    REQUIRE(negative.has_value());
    REQUIRE(negative.value().evaluate(Scalar{0}).velocity[0] == Scalar{-1});

    auto const leaving = create(Scalar{10}, Scalar{0}, one_ulp_above, 1);
    REQUIRE(leaving.has_value());
    REQUIRE(leaving.value().evaluate(leaving.value().duration()).velocity[0] == Scalar{1});

    REQUIRE(create(Scalar{10}, Scalar{1}, Scalar{0}, 0).has_value());
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

// --------------------------------------------------------------------------
// The out-parameter evaluation writes what the by-value one returns
// --------------------------------------------------------------------------
TEMPLATE_TEST_CASE("double_s: out-parameter evaluate is bit-identical to by-value evaluate",
                   "[traj][double_s]", float, double)
{
    using Scalar = TestType;
    using traj_t = ctrlpp::double_s_trajectory<Scalar>;

    SECTION("positive displacement, every segment")
    {
        require_bit_identical_over_horizon(built_profile<Scalar>(
            {.q0 = Scalar{0}, .q1 = Scalar{10}, .v_max = Scalar{5}, .a_max = Scalar{10}, .j_max = Scalar{100}}));
    }

    SECTION("negative displacement, where the sign frame is active")
    {
        auto const traj = built_profile<Scalar>(
            {.q0 = Scalar{10}, .q1 = Scalar{0}, .v_max = Scalar{5}, .a_max = Scalar{10}, .j_max = Scalar{100}});
        REQUIRE(traj.evaluate(traj.duration() / Scalar{2}).velocity[0] < Scalar{0});
        require_bit_identical_over_horizon(traj);
    }

    SECTION("negative displacement with a nonzero boundary velocity")
    {
        require_bit_identical_over_horizon(built_profile<Scalar>({.q0 = Scalar{0},
                                                                  .q1 = Scalar{-10},
                                                                  .v_max = Scalar{5},
                                                                  .a_max = Scalar{10},
                                                                  .j_max = Scalar{30},
                                                                  .v0 = Scalar{-1},
                                                                  .v1 = Scalar{0}}));
    }

    SECTION("each early return")
    {
        auto const standstill = built_profile<Scalar>(
            {.q0 = Scalar{5}, .q1 = Scalar{5}, .v_max = Scalar{5}, .a_max = Scalar{10}, .j_max = Scalar{100}});
        REQUIRE(standstill.duration() == Scalar{0});
        for (Scalar const t : {Scalar{-1}, Scalar{0}, Scalar{1}}) {
            require_bit_identical_at(standstill, t);
        }

        typename traj_t::config const cfg{.q0 = Scalar{1},
                                          .q1 = Scalar{-4},
                                          .v_max = Scalar{5},
                                          .a_max = Scalar{10},
                                          .j_max = Scalar{30},
                                          .v0 = Scalar{-1},
                                          .v1 = Scalar{-0.5}};
        auto const traj = built_profile<Scalar>(cfg);
        auto const T = traj.duration();
        for (Scalar const t : {T, Scalar{2} * T}) {
            Scalar q{}, dq{}, ddq{};
            traj.evaluate(t, q, dq, ddq);
            REQUIRE((q == cfg.q1 && dq == cfg.v1 && ddq == Scalar{0}));
            require_bit_identical_at(traj, t);
        }
        for (Scalar const t : {Scalar{0}, -T}) {
            Scalar q{}, dq{}, ddq{};
            traj.evaluate(t, q, dq, ddq);
            REQUIRE((q == cfg.q0 && dq == cfg.v0 && ddq == Scalar{0}));
            require_bit_identical_at(traj, t);
        }
    }
}

namespace
{

template <typename Scalar>
struct double_s_limits
{
    Scalar v_max;
    Scalar a_max;
    Scalar j_max;
};

template <typename Scalar>
struct sampled_excess
{
    std::int32_t acceleration;
    std::int32_t velocity;
};

template <typename Scalar>
void sample_into(ctrlpp::double_s_trajectory<Scalar> const& traj, double_s_limits<Scalar> const& l, Scalar t,
                 sampled_excess<Scalar>& excess)
{
    auto const point = traj.evaluate(t);
    excess.acceleration += (std::abs(point.acceleration[0]) > l.a_max) ? 1 : 0;
    excess.velocity += (std::abs(point.velocity[0]) > l.v_max) ? 1 : 0;
}

// The windows hold every representable instant within 64 steps of each cumulative switch, plus pred(T) and T.
template <typename Scalar>
auto excess_in_switch_windows(ctrlpp::double_s_trajectory<Scalar> const& traj, double_s_limits<Scalar> const& l)
    -> sampled_excess<Scalar>
{
    constexpr std::int32_t steps = 64;
    auto excess = sampled_excess<Scalar>{0, 0};
    auto center = Scalar{0};
    for (auto const segment : traj.phase_durations()) {
        center += segment;
        auto t = center;
        for (std::int32_t i = 0; i < steps; ++i) {
            t = std::nextafter(t, Scalar{0});
        }
        for (std::int32_t i = 0; i <= 2 * steps; ++i, t = std::nextafter(t, traj.duration())) {
            sample_into(traj, l, t, excess);
        }
    }
    sample_into(traj, l, std::nextafter(traj.duration(), Scalar{0}), excess);
    sample_into(traj, l, traj.duration(), excess);
    return excess;
}

template <typename Scalar>
auto excess_around_switches(ctrlpp::double_s_trajectory<Scalar> const& traj, double_s_limits<Scalar> const& l)
    -> sampled_excess<Scalar>
{
    constexpr auto grid_step = static_cast<Scalar>(1e-3);
    auto excess = excess_in_switch_windows(traj, l);
    for (std::int32_t i = 0; static_cast<Scalar>(i) * grid_step < traj.duration(); ++i) {
        sample_into(traj, l, static_cast<Scalar>(i) * grid_step, excess);
    }
    return excess;
}

template <typename Scalar>
void require_within_limits(double_s_limits<Scalar> const& l, Scalar q1, Scalar v0, Scalar v1)
{
    auto const created = ctrlpp::double_s_trajectory<Scalar>::create(
        {.q0 = 0, .q1 = q1, .v_max = l.v_max, .a_max = l.a_max, .j_max = l.j_max, .v0 = v0, .v1 = v1});
    REQUIRE(created.has_value());
    auto const excess = excess_around_switches(created.value(), l);
    CAPTURE(l.v_max, l.a_max, l.j_max, q1, v0, v1, excess.acceleration, excess.velocity);
    CHECK(excess.acceleration == 0);
    CHECK(excess.velocity == 0);
}

// j_max times the rounded quotient a_max / j_max rounds above a_max on this set.
template <typename Scalar>
auto round_up_limits() -> double_s_limits<Scalar>
{
    auto const a_max = std::is_same_v<Scalar, double> ? static_cast<Scalar>(0x1.fbfbfbdbfbfbfp-4) : static_cast<Scalar>(0.1);
    auto const j_max = std::is_same_v<Scalar, double> ? static_cast<Scalar>(1e6) : static_cast<Scalar>(3);
    return {Scalar{4} * a_max * a_max / j_max, a_max, j_max};
}

constexpr std::array<std::array<double, 2>, 8> sweep_rates{
    {{1, 1}, {10, 50}, {3, 7}, {1.3, 2.9}, {100, 1000}, {5, 50}, {20, 100}, {1e-2, 1e-1}}};

template <typename Scalar>
auto stepped_down(Scalar value, std::int32_t steps) -> Scalar
{
    for (std::int32_t i = 0; i < steps; ++i) {
        value = std::nextafter(value, Scalar{0});
    }
    return value;
}

template <typename Scalar>
auto stepped_up(Scalar value) -> Scalar
{
    return std::nextafter(value, std::numeric_limits<Scalar>::infinity());
}

}

TEMPLATE_TEST_CASE("Double-S: acceleration stays within the limit at every representable instant around each segment switch",
                   "[trajectory][double_s]", float, double)
{
    using Scalar = TestType;
    constexpr auto fast = double_s_limits<Scalar>{Scalar{40}, Scalar{20}, Scalar{100}};
    for (auto const sigma : {Scalar{1}, Scalar{-1}}) {
        require_within_limits(fast, sigma * Scalar{400}, Scalar{0}, sigma * fast.v_max);
        require_within_limits(fast, sigma * Scalar{400}, Scalar{0}, Scalar{0});
        require_within_limits(fast, sigma * Scalar{400}, sigma * fast.v_max, Scalar{0});
    }

    auto const round_up = round_up_limits<Scalar>();
    REQUIRE(round_up.j_max * (round_up.a_max / round_up.j_max) > round_up.a_max);
    auto const round_up_T_a = round_up.v_max / round_up.a_max + round_up.a_max / round_up.j_max;
    for (auto const sigma : {Scalar{1}, Scalar{-1}}) {
        require_within_limits(round_up, sigma * Scalar{10} * round_up.v_max * round_up_T_a, Scalar{0}, Scalar{0});
    }

    for (auto const& rates : sweep_rates) {
        auto const a_max = static_cast<Scalar>(rates[0]);
        auto const j_max = static_cast<Scalar>(rates[1]);
        auto v_max = a_max * a_max / j_max;
        for (std::int32_t i = 0; i < 8; ++i) {
            v_max = std::nextafter(v_max, Scalar{0});
        }
        for (std::int32_t k = -8; k <= 8; ++k, v_max = std::nextafter(v_max, std::numeric_limits<Scalar>::infinity())) {
            auto const T_a = v_max / a_max + a_max / j_max;
            require_within_limits(double_s_limits<Scalar>{v_max, a_max, j_max}, Scalar{10} * v_max * T_a, Scalar{0}, Scalar{0});
        }
    }
}

namespace
{

template <typename Scalar>
void require_switch_windows_within(ctrlpp::double_s_trajectory<Scalar> const& traj, double_s_limits<Scalar> const& l)
{
    auto const excess = excess_in_switch_windows(traj, l);
    CAPTURE(l.v_max, l.a_max, l.j_max, excess.acceleration, excess.velocity);
    CHECK(excess.acceleration == 0);
    CHECK(excess.velocity == 0);
}

// Full-mantissa magnitude over the given span of binary octaves, from raw engine words only.
template <typename Scalar>
auto seeded_magnitude(std::mt19937& engine, std::int32_t lowest_exponent, std::uint32_t octaves) -> Scalar
{
    constexpr auto fraction_bits = std::numeric_limits<Scalar>::digits - 1;
    auto const high = static_cast<std::uint64_t>(engine());
    auto const low = static_cast<std::uint64_t>(engine());
    auto const fraction = static_cast<Scalar>(((high << 32) | low) >> (64 - fraction_bits));
    auto const exponent = lowest_exponent + static_cast<std::int32_t>(engine() % octaves);
    return std::ldexp(Scalar{1} + std::ldexp(fraction, -fraction_bits), exponent);
}

template <typename Scalar>
auto seeded_fraction(std::mt19937& engine) -> Scalar
{
    return std::ldexp(static_cast<Scalar>(static_cast<std::int32_t>(engine() >> 8) - (1 << 23)), -23);
}

}

TEST_CASE("Double-S: a triangular ramp at its threshold holds its peak acceleration within the limit", "[trajectory][double_s]")
{
    auto const d = 0.50507005729686805;
    auto const l = double_s_limits<double>{2.0 * d, 0.029259168123834691, 0.0016950102405212067};
    auto const traj = built({.q0 = 0, .q1 = 1000, .v_max = l.v_max, .a_max = l.a_max, .j_max = l.j_max, .v0 = d, .v1 = 0});
    CHECK(traj.evaluate(traj.phase_durations()[0]).acceleration[0] <= l.a_max);
    require_switch_windows_within(traj, l);
}

// A jerk segment is sampled strictly before its stored end, so its largest sampled time is the predecessor of T_j.
TEMPLATE_TEST_CASE("Double-S: a velocity ramp near its triangular threshold never reaches above the acceleration limit",
                   "[trajectory][double_s]", float, double)
{
    using Scalar = TestType;
    std::mt19937 engine{std::is_same_v<Scalar, double> ? 501u : 502u};
    std::int32_t above = 0;
    auto first = std::array<Scalar, 3>{};
    for (std::int32_t pair = 0; pair < 4096; ++pair) {
        auto const a_max = seeded_magnitude<Scalar>(engine, -10, 20);
        auto const j_max = seeded_magnitude<Scalar>(engine, -10, 30);
        auto d = stepped_down(a_max * a_max / j_max, 16);
        for (std::int32_t k = -16; k <= 16; ++k, d = stepped_up(d)) {
            auto const ramp = ctrlpp::detail::make_velocity_ramp(d, a_max, j_max);
            if (!(ramp.a_lim <= a_max) || !(j_max * std::nextafter(ramp.T_j, Scalar{0}) <= a_max)) {
                first = (above == 0) ? std::array<Scalar, 3>{a_max, j_max, d} : first;
                ++above;
            }
        }
    }
    CAPTURE(first[0], first[1], first[2]);
    CHECK(above == 0);
}

TEMPLATE_TEST_CASE("Double-S: rest-to-rest profiles around the triangular threshold stay within the acceleration limit",
                   "[trajectory][double_s]", float, double)
{
    using Scalar = TestType;
    for (auto const& rates : sweep_rates) {
        auto const a_max = static_cast<Scalar>(rates[0]);
        auto const j_max = static_cast<Scalar>(rates[1]);
        auto v_max = stepped_down(a_max * a_max / j_max, 8);
        auto h = stepped_down(Scalar{2} * a_max * a_max * a_max / (j_max * j_max), 8);
        for (std::int32_t k = -8; k <= 8; ++k, v_max = stepped_up(v_max), h = stepped_up(h)) {
            auto const at_cruise_boundary = Scalar{2} * v_max * std::sqrt(v_max / j_max);
            auto const reaching = double_s_limits<Scalar>{Scalar{4} * a_max * a_max / j_max, a_max, j_max};
            for (auto const sigma : {Scalar{1}, Scalar{-1}}) {
                require_within_limits(double_s_limits<Scalar>{v_max, a_max, j_max}, sigma * at_cruise_boundary, Scalar{0}, Scalar{0});
                require_within_limits(reaching, sigma * h, Scalar{0}, Scalar{0});
            }
        }
    }
}

TEST_CASE("Double-S: the constant-deceleration segment of a long profile ends within the velocity limit", "[trajectory][double_s]")
{
    SECTION("float")
    {
        auto const v = static_cast<float>(0.0044597648084163666);
        auto const l = double_s_limits<float>{v, static_cast<float>(0.0015443010488525033), static_cast<float>(30.257682800292969)};
        auto const traj = built_profile<float>({.q0 = static_cast<float>(0.00033788720611482859),
                                                .q1 = static_cast<float>(-5.8095717430114746),
                                                .v_max = l.v_max,
                                                .a_max = l.a_max,
                                                .j_max = l.j_max,
                                                .v0 = static_cast<float>(-0.0035865840036422014),
                                                .v1 = v});
        CHECK(std::abs(traj.evaluate(std::nextafter(traj.duration(), 0.0f)).velocity[0]) <= v);
        require_switch_windows_within(traj, l);
    }

    SECTION("double")
    {
        auto const v = 0.086115804201940746;
        auto const l = double_s_limits<double>{v, 0.0037083513761102535, 4124533871.0118299};
        auto const traj = built({.q0 = 0, .q1 = 1015.7326884443336, .v_max = v, .a_max = l.a_max, .j_max = l.j_max, .v0 = 0, .v1 = -v});
        CHECK(std::abs(traj.evaluate(std::nextafter(traj.duration(), 0.0)).velocity[0]) <= v);
        require_switch_windows_within(traj, l);
    }
}

// Durations near 2^m resolve time to ulp(2^m), and each jerk segment lasts c of those units. The long-ramp variant
// stretches the acceleration ramp across a large part of the profile, which is where the constant-acceleration
// segment's switch is decided at a coarse resolution.
TEMPLATE_TEST_CASE("Double-S: long profiles with jerk segments near the time resolution stay within their limits",
                   "[trajectory][double_s]", float, double)
{
    using Scalar = TestType;
    constexpr std::array<std::array<double, 3>, 8> sets{{{1, 1, 1}, {5, 10, 50}, {2, 3, 7}, {0.7, 1.3, 2.9},
                                                         {10, 100, 1000}, {1, 5, 50}, {40, 20, 100}, {1e-3, 1e-2, 1e-1}}};
    constexpr auto digits = std::numeric_limits<Scalar>::digits;
    std::mt19937 engine{std::is_same_v<Scalar, double> ? 601u : 602u};
    auto failing = std::array<std::int32_t, 2>{0, 0};
    for (auto const& set : sets) {
        for (auto const m : {digits / 2, 3 * digits / 4}) {
            auto const span = std::ldexp(Scalar{1}, m);
            auto const resolution = stepped_up(span) - span;
            for (auto const long_ramp : {false, true}) {
                auto const v_max = static_cast<Scalar>(set[0]);
                auto const a_max = long_ramp ? Scalar{4} * v_max / span : static_cast<Scalar>(set[1]);
                for (auto const c : {Scalar{0.25}, Scalar{0.5}, Scalar{1}, Scalar{2}}) {
                    auto const l = double_s_limits<Scalar>{v_max, a_max, a_max / (c * resolution)};
                    for (auto const sigma : {Scalar{1}, Scalar{-1}}) {
                        for (auto const end : {Scalar{-1}, Scalar{1}}) {
                            for (std::int32_t k = 0; k < 8; ++k) {
                                auto const v0 = sigma * seeded_fraction<Scalar>(engine) * v_max;
                                auto const traj = built_profile<Scalar>({.q0 = Scalar{0},
                                                                         .q1 = sigma * v_max * span,
                                                                         .v_max = v_max,
                                                                         .a_max = a_max,
                                                                         .j_max = l.j_max,
                                                                         .v0 = v0,
                                                                         .v1 = end * sigma * v_max});
                                auto const excess = excess_in_switch_windows(traj, l);
                                failing[long_ramp ? 1 : 0] += (excess.velocity + excess.acceleration > 0) ? 1 : 0;
                            }
                        }
                    }
                }
            }
        }
    }
    CAPTURE(failing[0], failing[1]);
    CHECK(failing[0] == 0);
    CHECK(failing[1] == 0);
}
