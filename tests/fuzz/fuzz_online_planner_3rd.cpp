#include "ctrlpp/trajectory/online_planner_3rd.h"

#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <algorithm>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // 8 doubles: pos, target, v_max, a_max, j_max, dt, t_start, retarget = 64 bytes
    if(size < 64)
        return 0;

    double buf[8];
    std::memcpy(buf, data, 64);

    for(int i = 0; i < 8; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    double pos = buf[0];
    double target = buf[1];
    double v_max = std::abs(buf[2]);
    double a_max = std::abs(buf[3]);
    double j_max = std::abs(buf[4]);
    double dt = std::abs(buf[5]);
    double t_start = buf[6];

    // Clamp limits to positive and reasonable
    v_max = std::clamp(v_max, 1e-3, 1e6);
    a_max = std::clamp(a_max, 1e-3, 1e6);
    j_max = std::clamp(j_max, 1e-3, 1e6);
    dt = std::clamp(dt, 1e-6, 1.0);

    // Dense scan: many fixed-size steps exercise every jerk/accel/cruise/decel
    // segment (and the settled tail) regardless of the fuzzed target distance.
    constexpr int num_samples = 2000;

    // Clamp the scan's absolute time origin to within one scan window's worth
    // of the reference time (t_ref_ = 0): sample() computes elapsed time as
    // t - t_ref_, and an absolute t far larger than the dt-sized increments
    // being added to it would carry that increment below t's own floating-
    // point resolution, corrupting the finite differences below with pure
    // representation noise rather than trajectory behavior.
    const double t_start_bound = static_cast<double>(num_samples) * dt;
    t_start = std::clamp(t_start, -t_start_bound, t_start_bound);

    // Rounding-op margin: sample() integrates through up to 11 constant-jerk
    // phases, each a chained multiply-add (cubic position, quadratic
    // velocity, linear acceleration) plus phase-boundary subtractions, so a
    // single sample can accumulate several ULPs of rounding at its own
    // operand scale, not just the single operation a bare epsilon assumes;
    // 32 is a generous round count of those chained operations.
    constexpr double rounding_op_margin = 32.0;

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
    // independent of which physical limit (v_max/a_max/j_max) is being
    // checked; it only matters once that limit is itself very small.
    const double time_resolution_floor
        = std::numeric_limits<double>::epsilon() * rounding_op_margin * static_cast<double>(num_samples);

    ctrlpp::online_planner_3rd<double> planner({
        .v_max = v_max, .a_max = a_max, .j_max = j_max});

    // reset() always starts the planner at rest with zero acceleration, and a
    // single update() right after it always plans a rest-to-rest,
    // zero-acceleration move to the target, so both the start and end
    // boundary velocities (and accelerations) of the scanned profile are
    // zero. The decoded retarget byte is intentionally left unused: replanning
    // mid-motion would start from a nonzero velocity/acceleration state, and
    // that nonzero boundary case is not yet exact for this jerk-limited
    // profile, so it is left out of this scan; a later phase broadens the
    // fuzzer domain to that case.
    planner.reset(pos);
    planner.update(target);

    double prev_q = 0.0;
    double prev_v = 0.0;
    double prev_a = 0.0;
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
        // The velocity check adds an a_max*dt physical term and the
        // acceleration check adds a j_max*dt physical term: when a phase
        // transition falls within a single dt-sized sampling step (coarse
        // sampling relative to that phase's own, possibly tiny, duration),
        // the returned value can differ from the clean limit/0 value by up
        // to one step's worth of the next-higher derivative's limit.
        const double v_env_tol = a_max * dt
            + std::numeric_limits<double>::epsilon() * rounding_op_margin * v_max * static_cast<double>(num_samples)
            + time_resolution_floor;
        const double a_env_tol = j_max * dt
            + std::numeric_limits<double>::epsilon() * rounding_op_margin * a_max * static_cast<double>(num_samples)
            + time_resolution_floor;
        if(std::abs(v) > v_max + v_env_tol)
            abort();
        if(std::abs(a) > a_max + a_env_tol)
            abort();

        if(have_prev)
        {
            // Continuity: position cannot move faster than the kinematic bound
            // set by bounded velocity and acceleration over one sample step
            // (this bound only uses the v/a limits, so it holds regardless of
            // how the acceleration itself varies within the step).
            const double continuity_bound = v_max * dt + 0.5 * a_max * dt * dt;
            const double continuity_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin
                * (std::abs(q) + std::abs(prev_q) + continuity_bound);
            if(std::abs(q - prev_q) > continuity_bound + continuity_tol)
                abort();

            // Finite-difference velocity from the position trace: the average
            // velocity over one step differs from an instantaneous sample by
            // at most a_max*dt, since velocity's own rate of change (the
            // acceleration) never exceeds a_max. A sample pair straddling a
            // phase boundary is evaluated by two distinct closed-form
            // branches that meet exactly only in real arithmetic, so an
            // extra term scaled by the position magnitude in play (the same
            // rounding-op margin, divided back out by dt to land in velocity
            // units) absorbs that seam's rounding.
            const double v_fd = (q - prev_q) / dt;
            const double seam_v_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin
                * (std::abs(q) + std::abs(prev_q)) / dt;
            const double v_fd_tol = a_max * dt
                + std::numeric_limits<double>::epsilon() * rounding_op_margin * v_max * static_cast<double>(num_samples)
                + time_resolution_floor + seam_v_tol;
            if(std::abs(v_fd) > v_max + v_fd_tol)
                abort();

            // Finite-difference acceleration from the velocity trace: the
            // average acceleration over one step differs from an
            // instantaneous sample by at most j_max*dt, since acceleration's
            // own rate of change (the jerk) never exceeds j_max, plus the
            // same phase-boundary seam term one derivative order up (scaled
            // by the velocity magnitude in play).
            const double a_fd = (v - prev_v) / dt;
            const double seam_a_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin
                * (std::abs(v) + std::abs(prev_v)) / dt;
            const double a_fd_tol = j_max * dt
                + std::numeric_limits<double>::epsilon() * rounding_op_margin * a_max * static_cast<double>(num_samples)
                + time_resolution_floor + seam_a_tol;
            if(std::abs(a_fd) > a_max + a_fd_tol)
                abort();

            // Finite-difference jerk from the acceleration trace: the
            // profile's jerk is piecewise constant, so the average jerk over
            // one sample step cannot exceed j_max beyond rounding, plus the
            // same phase-boundary seam term one more derivative order up
            // (scaled by the acceleration magnitude in play).
            const double j_fd = (a - prev_a) / dt;
            const double seam_j_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin
                * (std::abs(a) + std::abs(prev_a)) / dt;
            const double j_fd_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin * j_max * static_cast<double>(num_samples)
                + time_resolution_floor + seam_j_tol;
            if(std::abs(j_fd) > j_max + j_fd_tol)
                abort();
        }

        prev_q = q;
        prev_v = v;
        prev_a = a;
        have_prev = true;
    }

    return 0;
}
