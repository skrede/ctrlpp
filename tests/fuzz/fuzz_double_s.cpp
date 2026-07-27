#include "ctrlpp/trajectory/double_s_trajectory.h"

#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <algorithm>

namespace
{

// Panels laid inside each of the seven phase segments for the swept-displacement
// quadrature. Simpson is exact on a velocity that is at most quadratic inside a
// segment, so once the panels stay inside one segment the quadrature carries no
// truncation error at all; a panel straddling a phase kink would carry one that
// reads as a violation on a correct profile.
constexpr int panels_per_segment = 4;

struct swept
{
    double integral = 0.0;
    double abs_area = 0.0;
    double back_time_exposure = 0.0;
    int panels = 0;
    bool usable = false;
};

// Kink-aligned integration of the reported velocity over the profile's own
// segments. This is the traversal oracle. A position sample is deliberately not
// used, at the ends of the interval or anywhere near them: the profile writes
// its final segment as an offset backwards from the commanded displacement, so
// it returns the target position by construction even when its velocity
// integrates to something else entirely, and a target built on that check could
// not find the defect class it exists to find.
swept integrate_velocity(const ctrlpp::double_s_trajectory<double>& traj)
{
    swept out;
    double t_start = 0.0;
    double widest_panel = 0.0;
    double last_length = 0.0;
    for(const auto& segment : traj.phase_durations())
    {
        if(!std::isfinite(segment) || segment < 0.0)
            return out;
        if(!(segment > 0.0))
            continue;

        const double dt = segment / static_cast<double>(panels_per_segment);
        widest_panel = std::max(widest_panel, dt);
        last_length = segment;
        for(int p = 0; p < panels_per_segment; ++p)
        {
            const double a = t_start + static_cast<double>(p) * dt;
            const double b = a + dt;
            const double m = a + 0.5 * dt;
            const double va = traj.evaluate(a).velocity(0);
            const double vm = traj.evaluate(m).velocity(0);
            const double vb = traj.evaluate(b).velocity(0);
            if(!std::isfinite(va) || !std::isfinite(vm) || !std::isfinite(vb))
                return out;
            const double panel = (dt / 6.0) * (va + 4.0 * vm + vb);
            if(!std::isfinite(panel))
                return out;
            out.integral += panel;
            out.abs_area += std::abs(panel);
            ++out.panels;
        }
        t_start += segment;
    }
    out.back_time_exposure = last_length + widest_panel;
    out.usable = std::isfinite(out.integral) && std::isfinite(out.abs_area) && out.panels > 0;
    return out;
}

}

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // Input layout, 8 little-endian IEEE-754 binary64 fields in this order,
    // minimum 64 bytes:
    //   q0, q1, v_max, a_max, j_max, v0, v1, stretch
    // The three limits are folded to their magnitude and clamped into a positive
    // range; the two boundary velocities are clamped into the range the decoded
    // velocity limit spans, so the fuzzer explores the whole boundary-velocity
    // family rather than a fixed slice of it; the stretch is folded to its
    // magnitude and clamped to at least one, and multiplies the constructed
    // duration to form the retiming request.
    if(size < 64)
        return 0;

    double buf[8];
    std::memcpy(buf, data, 64);

    for(int i = 0; i < 8; ++i)
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

    // Boundary velocities, in a range derived from the decoded velocity limit
    // rather than a fixed interval, so no bound of its own is introduced here.
    // The range is open at both ends: a boundary velocity that reaches the
    // velocity limit exactly leaves no admissible profile between the two, which
    // is outside the domain these profiles are defined on rather than a defect
    // in them. The neighbouring representable value is the largest speed still
    // inside that domain, so it is what the clamp uses.
    const double v_bound = std::nextafter(v_max, 0.0);
    const double v0 = std::clamp(buf[5], -v_bound, v_bound);
    const double v1 = std::clamp(buf[6], -v_bound, v_bound);

    // Retiming multiplier, on the same positive-range domain clamp the limits use.
    const double stretch = std::clamp(std::abs(buf[7]), 1.0, 1e6);

    // A command below the distance the fastest admissible transition between the
    // two boundary velocities already sweeps is not realizable by any
    // seven-segment shape, and the construction reports that rather than
    // returning a profile which does not traverse its own command. That is a
    // correct outcome, not a finding.
    auto created = ctrlpp::double_s_trajectory<double>::try_create(
        {.q0 = q0, .q1 = q1, .v_max = v_max, .a_max = a_max, .j_max = j_max, .v0 = v0, .v1 = v1});
    if(!created.has_value())
        return 0;
    ctrlpp::double_s_trajectory<double> traj = created.value();

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
    const double abs_h = std::abs(q1 - q0);
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
            // The deceleration segments form their position by subtracting the
            // distance still to run from the commanded displacement, so their
            // rounding lives at the scale of that displacement rather than at
            // the scale of the position value returned; when the move ends near
            // the origin those two differ by every digit.
            const double continuity_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin
                * (std::abs(q) + std::abs(prev_q) + abs_h + continuity_bound);
            if(std::abs(q - prev_q) > continuity_bound + continuity_tol)
                abort();

            // Finite-difference velocity from the position trace: the average
            // velocity over one step differs from an instantaneous sample by
            // at most a_max*dt, since velocity's own rate of change (the
            // acceleration) never exceeds a_max.
            const double v_fd = (q - prev_q) / dt;
            // The same displacement-scale rounding, divided back out by the step
            // to land in velocity units.
            const double seam_v_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin
                * (std::abs(q) + std::abs(prev_q) + abs_h) / dt;
            const double v_fd_tol = a_max * dt
                + std::numeric_limits<double>::epsilon() * rounding_op_margin * v_max * static_cast<double>(num_samples)
                + time_resolution_floor + seam_v_tol;
            if(std::abs(v_fd) > v_max + v_fd_tol)
                abort();

            // Finite-difference acceleration from the velocity trace: the
            // average acceleration over one step differs from an
            // instantaneous sample by at most j_max*dt, since acceleration's
            // own rate of change (the jerk) never exceeds j_max.
            const double a_fd = (v - prev_v) / dt;
            // A sample pair straddling a phase boundary is evaluated by two
            // distinct closed-form branches that meet exactly only in real
            // arithmetic, so one unit in the last place of the velocity
            // magnitude in play, divided back out by the step, has to be
            // absorbed here; at a large velocity over a short step that term
            // dwarfs the acceleration limit itself.
            const double seam_a_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin
                * (std::abs(v) + std::abs(prev_v)) / dt;
            const double a_fd_tol = j_max * dt
                + std::numeric_limits<double>::epsilon() * rounding_op_margin * a_max * static_cast<double>(num_samples)
                + time_resolution_floor + seam_a_tol;
            if(std::abs(a_fd) > a_max + a_fd_tol)
                abort();

            // Finite-difference jerk from the acceleration trace: the profile's
            // jerk is piecewise constant, so the average jerk over one sample
            // step cannot exceed j_max beyond rounding.
            const double j_fd = (a - prev_a) / dt;
            // The same phase-boundary seam term, one derivative order up.
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

    // -----------------------------------------------------------------------
    // Time-scaling leg
    // -----------------------------------------------------------------------
    if(!(T > 0.0))
        return 0;

    const double T_new = T * stretch;
    if(!std::isfinite(T_new) || !(T_new > T))
        return 0;

    ctrlpp::double_s_trajectory<double> retimed = traj;
    const auto rescaled = retimed.rescale_to(T_new);
    if(!rescaled.has_value())
    {
        // A typed rejection is a correct outcome for a duration no admissible
        // scale realizes, so it is not a finding. What would be a finding is a
        // rejection that left the profile changed anyway.
        if(retimed.duration() != T)
            abort();
        return 0;
    }

    const double T_scaled = retimed.duration();
    if(!std::isfinite(T_scaled) || !(T_scaled > 0.0))
        abort();

    // Phase bookkeeping: no phase may be negative, and the reported duration
    // must be what the phases add up to.
    double segment_sum = 0.0;
    for(const auto& segment : retimed.phase_durations())
    {
        if(!std::isfinite(segment) || segment < 0.0)
            abort();
        segment_sum += segment;
    }
    if(std::abs(segment_sum - T_scaled)
       > std::numeric_limits<double>::epsilon() * rounding_op_margin * T_scaled)
        abort();

    // Swept displacement, by kink-aligned quadrature. Two terms: the
    // quadrature's own rounding floor at the scale of the absolute area its
    // panels accumulated, and the deceleration branch's backward-time recovery,
    // which forms its local time by subtracting from the duration and so carries
    // up to one unit in the last place of that duration, converted into velocity
    // through the slew there. The weight exposed to it is the final segment plus
    // one panel of its predecessor, whose right endpoint is the shared boundary.
    const swept area = integrate_velocity(retimed);
    if(!area.usable)
        return 0;

    const double h_signed = q1 - q0;
    const double area_scale = std::max(area.abs_area, std::abs(h_signed));
    const double displacement_tol
        = std::numeric_limits<double>::epsilon() * rounding_op_margin
        * (static_cast<double>(area.panels) * area_scale + a_max * T_scaled * area.back_time_exposure);
    if(!std::isfinite(displacement_tol))
        return 0;
    if(std::abs(area.integral - h_signed) > displacement_tol)
        abort();

    // Boundary velocities, a step inside each end. The step is the profile's own
    // first (respectively last) nonempty segment, and the excursion across it
    // cannot exceed the acceleration limit times its length whatever segment
    // shape it has.
    double first_length = 0.0;
    double last_length = 0.0;
    for(const auto& segment : retimed.phase_durations())
    {
        if(!(segment > 0.0))
            continue;
        if(!(first_length > 0.0))
            first_length = segment;
        last_length = segment;
    }

    const double v_ceiling = std::max({v_max, std::abs(v0), std::abs(v1)});
    const double v_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin
        * (v_ceiling + a_max * T_scaled);
    if(!std::isfinite(v_tol))
        return 0;

    if(first_length > 0.0)
    {
        const double v_start = retimed.evaluate(first_length).velocity(0);
        if(!std::isfinite(v_start))
            return 0;
        if(std::abs(v_start - v0) > a_max * first_length + v_tol)
            abort();
    }
    if(last_length > 0.0)
    {
        const double v_end = retimed.evaluate(T_scaled - last_length).velocity(0);
        if(!std::isfinite(v_end))
            return 0;
        if(std::abs(v_end - v1) > a_max * last_length + v_tol)
            abort();
    }

    // Kinematic limits over the retimed profile, on a grid laid over the
    // requested duration rather than the realized one. The retiming divides each
    // limit by a power of the scale it solves for, so the original limits stay
    // valid envelopes for the retimed profile.
    constexpr int scaling_samples = 64;
    const double scaled_dt = T_new / static_cast<double>(scaling_samples);
    for(int i = 0; i <= scaling_samples; ++i)
    {
        const auto point = retimed.evaluate(static_cast<double>(i) * scaled_dt);
        const double v_sample = point.velocity(0);
        const double a_sample = point.acceleration(0);
        if(!std::isfinite(v_sample) || !std::isfinite(a_sample))
            return 0;
        // The same backward-time margin the boundary check carries: near the end
        // of a long move one unit in the last place of the duration converts,
        // through the slew there, into a velocity uncertainty that can exceed
        // the velocity limit itself.
        if(std::abs(v_sample) > v_ceiling + v_tol)
            abort();
        if(std::abs(a_sample)
           > a_max + std::numeric_limits<double>::epsilon() * rounding_op_margin * a_max)
            abort();
    }

    // Realized duration. The scale is either a single quotient of the two
    // durations, when the axis starts and ends at rest, or the endpoint of a
    // bracket halved until its midpoint lands on an endpoint. Either way it is
    // settled to within one unit in the last place, and within a fixed segment
    // shape the duration is proportional to the reciprocal of the scale, so the
    // duration inherits that unit at its own magnitude; the remaining budget
    // covers the scaled limits, the phase durations derived from them, and their
    // sum.
    const double duration_tol
        = std::numeric_limits<double>::epsilon() * rounding_op_margin * T_new;
    if(std::isfinite(duration_tol) && std::abs(T_scaled - T_new) > duration_tol)
        abort();

    return 0;
}
