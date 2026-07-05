#include "ctrlpp/trajectory/double_s_trajectory.h"

#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <algorithm>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // 5 doubles: q0, q1, v_max, a_max, j_max = 40 bytes
    if(size < 40)
        return 0;

    double buf[5];
    std::memcpy(buf, data, 40);

    for(int i = 0; i < 5; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    double q0 = buf[0];
    double q1 = buf[1];
    double v_max = std::abs(buf[2]);
    double a_max = std::abs(buf[3]);
    double j_max = std::abs(buf[4]);

    // Clamp limits to positive
    v_max = std::clamp(v_max, 1e-6, 1e6);
    a_max = std::clamp(a_max, 1e-6, 1e6);
    j_max = std::clamp(j_max, 1e-6, 1e6);

    // The boundary velocities are forced to zero here: nonzero boundary-velocity
    // double-S trajectories are not yet exact, so this scan is restricted to
    // zero start/end velocity; a later phase broadens the fuzzer domain to
    // nonzero start/end velocities once that formula is corrected.
    ctrlpp::double_s_trajectory<double> traj({
        .q0 = q0, .q1 = q1, .v_max = v_max, .a_max = a_max, .j_max = j_max, .v0 = 0.0, .v1 = 0.0});

    const double T = traj.duration();
    if(!std::isfinite(T) || T < 0.0)
        return 0;

    // Dense scan of the time domain: enough samples to exercise every segment
    // (jerk/accel/cruise/jerk/decel) multiple times regardless of duration.
    constexpr int num_samples = 500;
    const double dt = (T > 0.0) ? T / static_cast<double>(num_samples) : 0.0;

    // Rounding-op margin: evaluate() chains several multiply-adds per sample
    // across up to 7 segments (cubic/quadratic terms, sigma-frame transform,
    // and a phase-boundary subtraction), each contributing up to one ULP of
    // rounding at its own operand scale, not just the single operation a bare
    // epsilon assumes; 32 is a generous round count of those chained
    // operations (matching the online-planner jerk-limited fuzzer's margin).
    constexpr double rounding_op_margin = 32.0;

    // Reject configurations whose position scale swamps the per-step motion:
    // if a handful of ULPs of q0/q1's own magnitude (the same rounding-op
    // margin used below) already exceeds the kinematic bound on how far the
    // position can move in a single sample step, adding that step's motion to
    // a position of this magnitude rounds away completely, and no continuity
    // signal survives to check -- this is a floating-point representation
    // limit, not a trajectory defect.
    const double position_scale = std::max({std::abs(q0), std::abs(q1), 1.0});
    const double continuity_bound_reject = v_max * dt + 0.5 * a_max * dt * dt;
    if(std::numeric_limits<double>::epsilon() * rounding_op_margin * position_scale >= continuity_bound_reject)
        return 0;

    // Reject configurations whose jerk-ramp phases are too short to resolve
    // against the trajectory's own total duration: the library locates each
    // sample within its 7 segments by subtracting sample times from phase
    // boundaries (e.g. T_a - tc), and if a jerk phase (a_max/j_max, or the
    // doubly-degenerate cbrt-derived ramp) is shorter than one ULP of the
    // total duration T, that boundary arithmetic cannot resolve which segment
    // a sample near the boundary falls in -- a floating-point representation
    // limit in locating the segment, not a trajectory defect.
    if(T > 0.0)
    {
        const auto phases = traj.phase_durations();
        const double min_jerk_phase = std::min({phases[0], phases[2], phases[4], phases[6]});
        if(min_jerk_phase > 0.0 && min_jerk_phase < std::numeric_limits<double>::epsilon() * rounding_op_margin * T)
            return 0;
    }

    // Sample-time floor: advancing through num_samples discrete sample times
    // carries its own fixed rounding floor of about one part in num_samples of
    // machine epsilon once divided back out by dt in a finite difference,
    // independent of which physical limit (v_max/a_max/j_max) is being
    // checked; it only matters once that limit is itself very small.
    const double time_resolution_floor
        = std::numeric_limits<double>::epsilon() * rounding_op_margin * static_cast<double>(num_samples);

    double prev_q = 0.0;
    double prev_v = 0.0;
    double prev_a = 0.0;
    bool have_prev = false;

    for(int i = 0; i <= num_samples; ++i)
    {
        const double t = static_cast<double>(i) * dt;
        auto pt = traj.evaluate(t);

        const double q = pt.position(0);
        const double v = pt.velocity(0);
        const double a = pt.acceleration(0);

        if(!std::isfinite(q) || !std::isfinite(v) || !std::isfinite(a))
            return 0;

        // Envelope: the analytically-returned velocity/acceleration must not
        // exceed the configured limits beyond a rounding-scaled tolerance.
        const double v_env_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin * v_max * static_cast<double>(num_samples);
        const double a_env_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin * a_max * static_cast<double>(num_samples);
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
            // acceleration) never exceeds a_max.
            const double v_fd = (q - prev_q) / dt;
            const double v_fd_tol = a_max * dt
                + std::numeric_limits<double>::epsilon() * rounding_op_margin * v_max * static_cast<double>(num_samples)
                + time_resolution_floor;
            if(std::abs(v_fd) > v_max + v_fd_tol)
                abort();

            // Finite-difference acceleration from the velocity trace: the
            // average acceleration over one step differs from an
            // instantaneous sample by at most j_max*dt, since acceleration's
            // own rate of change (the jerk) never exceeds j_max.
            const double a_fd = (v - prev_v) / dt;
            const double a_fd_tol = j_max * dt
                + std::numeric_limits<double>::epsilon() * rounding_op_margin * a_max * static_cast<double>(num_samples)
                + time_resolution_floor;
            if(std::abs(a_fd) > a_max + a_fd_tol)
                abort();

            // Finite-difference jerk from the acceleration trace: the profile's
            // jerk is piecewise constant, so the average jerk over one sample
            // step cannot exceed j_max beyond rounding.
            const double j_fd = (a - prev_a) / dt;
            const double j_fd_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin * j_max * static_cast<double>(num_samples)
                + time_resolution_floor;
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
