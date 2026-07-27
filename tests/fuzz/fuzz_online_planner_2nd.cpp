#include "ctrlpp/trajectory/online_planner_2nd.h"

#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <algorithm>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // 7 doubles: pos, vel, target, v_max, a_max, dt, t_start = 56 bytes
    if(size < 56)
        return 0;

    double buf[7];
    std::memcpy(buf, data, 56);

    for(int i = 0; i < 7; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    double pos = buf[0];
    double target = buf[2];
    double v_max = std::abs(buf[3]);
    double a_max = std::abs(buf[4]);
    double dt = std::abs(buf[5]);
    double t_start = buf[6];

    // Clamp limits to positive and reasonable
    v_max = std::clamp(v_max, 1e-3, 1e6);
    a_max = std::clamp(a_max, 1e-3, 1e6);
    dt = std::clamp(dt, 1e-6, 1.0);

    // Dense scan: many fixed-size steps exercise the accel/cruise/decel
    // phases (and the settled tail) regardless of the fuzzed target distance.
    constexpr int num_samples = 2000;

    // Clamp the scan's absolute time origin to within one scan window's worth
    // of the reference time (t_ref_ = 0): sample() computes elapsed time as
    // t - t_ref_, and an absolute t far larger than the dt-sized increments
    // being added to it would carry that increment below t's own floating-
    // point resolution, corrupting the finite differences below with pure
    // representation noise rather than trajectory behavior.
    const double t_start_bound = static_cast<double>(num_samples) * dt;
    t_start = std::clamp(t_start, -t_start_bound, t_start_bound);

    // Rounding-op margin: each sample() call chains several multiply-adds and
    // a phase-boundary subtraction (backward time in the deceleration phase),
    // each contributing up to one ULP of rounding at its own operand scale,
    // not just the single operation a bare epsilon assumes; 16 is a generous
    // round count of those chained operations.
    constexpr double rounding_op_margin = 16.0;

    // Reject configurations whose position scale swamps the per-step motion:
    // if a handful of ULPs of pos/target's own magnitude (the same
    // rounding-op margin used below) already exceeds the kinematic bound on
    // how far the position can move in a single sample step, adding that
    // step's motion to a position of this magnitude rounds away completely,
    // and no continuity signal survives to check -- this is a floating-point
    // representation limit, not a planner defect.
    const double position_scale = std::max({std::abs(pos), std::abs(target), 1.0});
    const double continuity_bound_reject = v_max * dt + 0.5 * a_max * dt * dt;
    if(std::numeric_limits<double>::epsilon() * rounding_op_margin * position_scale >= continuity_bound_reject)
        return 0;

    // Sample-time floor: advancing through num_samples discrete sample times
    // carries its own fixed rounding floor of about one part in num_samples of
    // machine epsilon once divided back out by dt in a finite difference,
    // independent of which physical limit (v_max/a_max) is being checked; it
    // only matters once that limit is itself very small.
    const double time_resolution_floor
        = std::numeric_limits<double>::epsilon() * rounding_op_margin * static_cast<double>(num_samples);

    // Construction goes through create: a limit outside the finite,
    // strictly positive domain is rejected by design (each limit divides in
    // the planner math), and that rejection is the correct behavior for such
    // an input, not a crash. The clamps above keep the fuzzed limits inside
    // the domain, so the scan below runs for every input that reaches here.
    auto planner_result
        = ctrlpp::online_planner_2nd<double>::create({.v_max = v_max, .a_max = a_max});
    if(!planner_result.has_value())
        return 0;
    auto& planner = *planner_result;

    // reset() always starts the planner at rest, and a single update() right
    // after it always plans a rest-to-rest move to the target, so both the
    // start and end boundary velocities of the scanned profile are zero (the
    // decoded velocity byte is intentionally left unused). Nonzero
    // boundary-velocity replanning -- e.g. retargeting while still moving --
    // is not yet exact for this trapezoidal-family profile, so it is left out
    // of this scan; a later phase broadens the fuzzer domain to that case.
    planner.reset(pos);
    planner.update(target);

    double prev_q = 0.0;
    double prev_v = 0.0;
    bool have_prev = false;
    double t = t_start;

    for(int i = 0; i < num_samples; ++i)
    {
        t += dt;
        auto pt = planner.sample(t);

        const double q = pt.position(0);
        const double v = pt.velocity(0);
        const double a = pt.acceleration(0);

        if(!std::isfinite(q) || !std::isfinite(v) || !std::isfinite(a))
            return 0;

        // Envelope: the analytically-returned velocity/acceleration must not
        // exceed the configured limits beyond a rounding-scaled tolerance.
        // The velocity check adds an a_max*dt physical term: when the
        // acceleration-to-cruise phase transition falls within a single
        // dt-sized sampling step (coarse sampling relative to that phase's
        // own, possibly tiny, duration), the returned v can differ from the
        // clean v_max/0 value by up to one step's worth of acceleration.
        const double v_env_tol = a_max * dt
            + std::numeric_limits<double>::epsilon() * rounding_op_margin * v_max * static_cast<double>(num_samples)
            + time_resolution_floor;
        const double a_env_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin * a_max * static_cast<double>(num_samples)
            + time_resolution_floor;
        if(std::abs(v) > v_max + v_env_tol)
            abort();
        if(std::abs(a) > a_max + a_env_tol)
            abort();

        if(have_prev)
        {
            // Continuity: position cannot move faster than the kinematic bound
            // set by bounded velocity and acceleration over one sample step.
            const double continuity_bound = v_max * dt + 0.5 * a_max * dt * dt;
            const double continuity_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin
                * (std::abs(q) + std::abs(prev_q) + continuity_bound);
            if(std::abs(q - prev_q) > continuity_bound + continuity_tol)
                abort();

            // Finite-difference velocity from the position trace: within a
            // single accel/cruise/decel phase, position is quadratic in t so
            // the forward difference matches the true derivative up to a
            // term bounded by one acceleration-phase jump (at most a_max)
            // times the step; a small eps floor absorbs rounding. A sample
            // pair straddling a phase boundary is evaluated by two distinct
            // closed-form branches that meet exactly only in real arithmetic,
            // so an extra term scaled by the position magnitude in play (the
            // same rounding-op margin, divided back out by dt to land in
            // velocity units) absorbs that seam's rounding.
            const double v_fd = (q - prev_q) / dt;
            const double seam_v_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin
                * (std::abs(q) + std::abs(prev_q)) / dt;
            const double v_fd_tol = a_max * dt
                + std::numeric_limits<double>::epsilon() * rounding_op_margin * v_max * static_cast<double>(num_samples)
                + time_resolution_floor + seam_v_tol;
            if(std::abs(v_fd) > v_max + v_fd_tol)
                abort();

            // Finite-difference acceleration from the returned velocity trace:
            // velocity is piecewise linear in t (constant acceleration per
            // phase), so this difference is a convex combination of at most
            // two phases' accelerations and is bounded by a_max up to
            // rounding, plus the same phase-boundary seam term one derivative
            // order up (scaled by the velocity magnitude in play).
            const double a_fd = (v - prev_v) / dt;
            const double seam_a_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin
                * (std::abs(v) + std::abs(prev_v)) / dt;
            const double a_fd_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin * a_max * static_cast<double>(num_samples)
                + time_resolution_floor + seam_a_tol;
            if(std::abs(a_fd) > a_max + a_fd_tol)
                abort();
        }

        prev_q = q;
        prev_v = v;
        have_prev = true;
    }

    return 0;
}
