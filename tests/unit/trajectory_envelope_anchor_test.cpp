// This anchor dense-scans the trapezoidal and double-S velocity profiles
// over randomized configurations that include nonzero initial and final
// velocities -- the one part of each profile's public configuration surface
// that the equivalent fuzz oracles do not yet exercise (those are
// deliberately restricted to zero boundary velocity until this behavior is
// corrected).
//
// Two properties are checked purely from each profile's own reported
// position/velocity trace, without reference to how that trace was
// produced internally:
//
//   - Continuity: the position cannot move, between two nearby sample
//     times, faster than the kinematic bound set by the configured
//     velocity limit (plus the extra distance one sample step of the
//     configured acceleration limit could contribute). Any larger jump is a
//     genuine discontinuity in the reported position trace, not a sampling
//     artifact.
//   - Envelope via finite difference: the reported velocity itself cannot
//     change, between two nearby sample times, faster than the configured
//     acceleration limit; and for the jerk-limited profile, the reported
//     acceleration cannot change faster than the configured jerk limit.
//
// The trapezoidal profile currently violates the continuity property once
// boundary velocities are nonzero: the two phases straddling the transition
// out of the cruise phase are stitched from inconsistent position
// expressions, producing an interior jump even though the reported velocity
// itself remains continuous through that same transition. The double-S
// profile currently violates the finite-difference acceleration envelope at
// its own boundaries: it silently computes its entire phase timing as if
// both boundary velocities were zero, so the reported velocity leaps from
// the configured nonzero boundary value to the zero-boundary profile's own
// value within a single sample step. Both sections are held green with
// `[!shouldfail]` until their respective profiles are corrected to handle
// nonzero boundary velocities exactly.

#include "ctrlpp/trajectory/double_s_trajectory.h"
#include "ctrlpp/trajectory/trapezoidal_trajectory.h"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <optional>

using namespace ctrlpp;

namespace
{

// Rounding-op margin: each evaluate() call chains several multiply-adds and,
// for the profiles compared here, a phase-boundary subtraction (backward
// time from the end of the move); each contributes up to one ULP of
// rounding at its own operand scale, not just a single bare epsilon. This
// margin matches the one used by this suite's equivalent fuzz oracle for
// the same closed-form profile arithmetic.
constexpr double rounding_op_margin = 16.0;

/// Dense-scan continuity and finite-difference velocity/acceleration/jerk
/// envelope check, shared by both profile types. `j_max` is empty for
/// profiles (trapezoidal) that do not bound jerk.
template <typename Trajectory>
void check_continuity_and_envelope(const Trajectory& trajectory, double v_max, double a_max, std::optional<double> j_max)
{
    const double T = trajectory.duration();
    if(!(T > 0.0))
        return; // degenerate zero-displacement profile: nothing to scan

    constexpr int num_samples = 400;
    const double dt = T / static_cast<double>(num_samples);
    const double eps = std::numeric_limits<double>::epsilon();

    auto p0 = trajectory.evaluate(0.0);
    double prev_q = p0.position(0);
    double prev_v = p0.velocity(0);
    double prev_a = p0.acceleration(0);

    for(int i = 1; i <= num_samples; ++i)
    {
        const double t = static_cast<double>(i) * dt;
        const auto pt = trajectory.evaluate(t);
        const double q = pt.position(0);
        const double v = pt.velocity(0);
        const double a = pt.acceleration(0);

        REQUIRE(std::isfinite(q));
        REQUIRE(std::isfinite(v));
        REQUIRE(std::isfinite(a));

        const double continuity_bound = v_max * dt + 0.5 * a_max * dt * dt;
        const double continuity_tol = eps * rounding_op_margin * (std::abs(q) + std::abs(prev_q) + continuity_bound);
        CAPTURE(t, q, prev_q, continuity_bound);
        REQUIRE(std::abs(q - prev_q) <= continuity_bound + continuity_tol);

        const double v_fd = std::abs(v - prev_v) / dt;
        const double a_env_tol = eps * rounding_op_margin * a_max * static_cast<double>(num_samples);
        CAPTURE(t, v, prev_v, v_fd);
        REQUIRE(v_fd <= a_max + a_env_tol);

        if(j_max.has_value())
        {
            const double a_fd = std::abs(a - prev_a) / dt;
            const double j_env_tol = eps * rounding_op_margin * (*j_max) * static_cast<double>(num_samples);
            CAPTURE(t, a, prev_a, a_fd);
            REQUIRE(a_fd <= *j_max + j_env_tol);
        }

        prev_q = q;
        prev_v = v;
        prev_a = a;
    }
}

struct random_config
{
    double q0, q1, v_max, a_max, j_max, v0, v1;
};

/// Six randomized configurations, fixed-seeded for reproducibility, whose
/// boundary velocities are a random nonzero fraction of v_max in the
/// direction of the move.
auto make_random_configs() -> std::array<random_config, 6>
{
    std::mt19937 gen(20260705);
    std::uniform_real_distribution<double> q_dist(-20.0, 20.0);
    std::uniform_real_distribution<double> v_max_dist(1.0, 5.0);
    std::uniform_real_distribution<double> a_max_dist(0.5, 3.0);
    std::uniform_real_distribution<double> j_max_dist(1.0, 8.0);
    std::uniform_real_distribution<double> boundary_frac_dist(0.05, 0.6);

    std::array<random_config, 6> configs{};
    for(auto& cfg : configs)
    {
        double q0 = q_dist(gen);
        double q1 = q_dist(gen);
        while(std::abs(q1 - q0) < 1.0)
            q1 = q_dist(gen);
        const double v_max = v_max_dist(gen);
        const double a_max = a_max_dist(gen);
        const double j_max = j_max_dist(gen);
        const double sign = (q1 > q0) ? 1.0 : -1.0;
        const double v0 = sign * boundary_frac_dist(gen) * v_max;
        const double v1 = sign * boundary_frac_dist(gen) * v_max;
        cfg = {q0, q1, v_max, a_max, j_max, v0, v1};
    }
    return configs;
}

} // namespace

TEST_CASE("trapezoidal profile continuity and velocity/acceleration envelope hold over randomized nonzero-boundary-velocity configurations",
    "[trajectory][anchor]")
{
    for(const auto& cfg : make_random_configs())
    {
        CAPTURE(cfg.q0, cfg.q1, cfg.v_max, cfg.a_max, cfg.v0, cfg.v1);
        trapezoidal_trajectory<double> trajectory({
            .q0 = cfg.q0, .q1 = cfg.q1, .v_max = cfg.v_max, .a_max = cfg.a_max, .v0 = cfg.v0, .v1 = cfg.v1});
        check_continuity_and_envelope(trajectory, cfg.v_max, cfg.a_max, std::nullopt);
    }
}

TEST_CASE("double-S profile continuity and velocity/acceleration/jerk envelope hold over randomized nonzero-boundary-velocity configurations",
    "[trajectory][anchor][!shouldfail]")
{
    for(const auto& cfg : make_random_configs())
    {
        CAPTURE(cfg.q0, cfg.q1, cfg.v_max, cfg.a_max, cfg.j_max, cfg.v0, cfg.v1);
        double_s_trajectory<double> trajectory({
            .q0 = cfg.q0, .q1 = cfg.q1, .v_max = cfg.v_max, .a_max = cfg.a_max, .j_max = cfg.j_max, .v0 = cfg.v0, .v1 = cfg.v1});
        check_continuity_and_envelope(trajectory, cfg.v_max, cfg.a_max, cfg.j_max);
    }
}
