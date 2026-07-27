#include "ctrlpp/trajectory/trapezoidal_trajectory.h"

#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <algorithm>

namespace
{

// Panels laid inside each phase segment for the swept-displacement quadrature.
// Simpson integrates the piecewise-linear velocity exactly once the panels stay
// inside one segment, so the count buys no accuracy; a panel straddling a phase
// kink would carry a truncation error that reads as a violation on a correct
// profile.
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
swept integrate_velocity(const ctrlpp::trapezoidal_trajectory<double>& traj)
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

// Amplification the cruise-velocity solve applies to its own roundings. Each of
// the three shapes reaches its root through a difference of two nearly equal
// quantities, and the ratio of the operands entering that difference to the
// difference itself is what the realized duration inherits.
double solve_conditioning(double v_cruise, double a, double h, double v0, double v1, double T_new)
{
    const double v_lo = std::min(v0, v1);
    const double v_hi = std::max(v0, v1);
    const double v_sum_sq = v0 * v0 + v1 * v1;

    const auto ratio = [](double operands, double residual) {
        const double scale = std::abs(residual);
        return (scale > 0.0) ? std::abs(operands) / scale : 0.0;
    };

    if(v_cruise > v_lo && v_cruise < v_hi)
    {
        const double ramp_distance = (v_hi - v_lo) * (v_hi + v_lo) / (2.0 * a);
        return std::max(1.0, ratio(std::max(h, v_sum_sq / (2.0 * a)), h - ramp_distance));
    }

    const bool plateau = (v_cruise >= v_hi);
    const double b = plateau ? ((v0 + v1) + a * T_new) : (a * T_new - (v0 + v1));
    const double c = plateau ? (a * h + v_sum_sq / 2.0) : (v_sum_sq / 2.0 - a * h);
    return std::max({1.0, ratio(a * h + v_sum_sq / 2.0, c), ratio(b * b + 4.0 * std::abs(c), b * b - 4.0 * c)});
}

}

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // Input layout, 7 little-endian IEEE-754 binary64 fields in this order,
    // minimum 56 bytes:
    //   q0, q1, v_max, a_max, v0, v1, stretch
    // The two limits are folded to their magnitude and clamped into a positive
    // range; the two boundary velocities are clamped into the range the decoded
    // velocity limit spans, so the fuzzer explores the whole boundary-velocity
    // family rather than a fixed slice of it; the stretch is folded to its
    // magnitude and clamped to at least one, and multiplies the constructed
    // duration to form the retiming request.
    if(size < 56)
        return 0;

    double buf[7];
    std::memcpy(buf, data, 56);

    for(int i = 0; i < 7; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    double q0 = buf[0];
    double q1 = buf[1];
    double v_max = std::abs(buf[2]);
    double a_max = std::abs(buf[3]);

    // Clamp limits to positive
    v_max = std::clamp(v_max, 1e-6, 1e6);
    a_max = std::clamp(a_max, 1e-6, 1e6);

    // Boundary velocities, in a range derived from the decoded velocity limit
    // rather than a fixed interval, so no bound of its own is introduced here.
    // The range is open at both ends: a boundary velocity that reaches the
    // velocity limit exactly leaves no admissible profile between the two, which
    // is outside the domain these profiles are defined on rather than a defect
    // in them. The neighbouring representable value is the largest speed still
    // inside that domain, so it is what the clamp uses.
    const double v_bound = std::nextafter(v_max, 0.0);
    const double v0 = std::clamp(buf[4], -v_bound, v_bound);
    const double v1 = std::clamp(buf[5], -v_bound, v_bound);

    // Retiming multiplier, on the same positive-range domain clamp the limits use.
    const double stretch = std::clamp(std::abs(buf[6]), 1.0, 1e6);

    ctrlpp::trapezoidal_trajectory<double> traj({
        .q0 = q0, .q1 = q1, .v_max = v_max, .a_max = a_max, .v0 = v0, .v1 = v1});

    const double T = traj.duration();
    if(!std::isfinite(T) || T < 0.0)
        return 0;

    // Effective acceleration magnitude. When the commanded displacement is too
    // short for the two boundary velocities at the commanded acceleration, the
    // construction raises that acceleration to the value which makes the two
    // ramps exactly cover the displacement, so the envelope the trace has to
    // respect is that raised value, not the decoded one.
    const double abs_h = std::abs(q1 - q0);
    const double ramp_only = std::abs(v0 * v0 - v1 * v1) / 2.0;
    const double a_eff = (abs_h > 0.0 && a_max * abs_h < ramp_only)
        ? (ramp_only / abs_h + std::numeric_limits<double>::epsilon())
        : a_max;
    if(!std::isfinite(a_eff))
        return 0;

    // Dense scan of the time domain: enough samples to exercise every phase
    // (accel/cruise/decel) multiple times regardless of the fuzzed duration.
    constexpr int num_samples = 500;
    const double dt = (T > 0.0) ? T / static_cast<double>(num_samples) : 0.0;

    // Rounding-op margin: each evaluate() call chains several multiply-adds
    // and a phase-boundary subtraction (backward time in the deceleration
    // phase), each contributing up to one ULP of rounding at its own operand
    // scale, not just the single operation a bare epsilon assumes; 16 is a
    // generous round count of those chained operations (matching the
    // second-order online planner's margin).
    constexpr double rounding_op_margin = 16.0;

    // Reject configurations whose position scale swamps the per-step motion:
    // if a handful of ULPs of q0/q1's own magnitude (the same rounding-op
    // margin used below) already exceeds the kinematic bound on how far the
    // position can move in a single sample step, adding that step's motion to
    // a position of this magnitude rounds away completely, and no continuity
    // signal survives to check -- this is a floating-point representation
    // limit, not a trajectory defect.
    const double position_scale = std::max({std::abs(q0), std::abs(q1), 1.0});
    const double continuity_bound_reject = v_max * dt + 0.5 * a_eff * dt * dt;
    if(std::numeric_limits<double>::epsilon() * rounding_op_margin * position_scale >= continuity_bound_reject)
        return 0;

    // No admissible speed exceeds the larger of the velocity limit and the two
    // commanded boundary velocities.
    const double v_ceiling = std::max({v_max, std::abs(v0), std::abs(v1)});

    // Conditioning of a phase duration, which the seam between two phases
    // inherits. The cruise velocity is a difference of two nearly equal
    // quantities against the boundary velocities, so its absolute error lives at
    // the velocity scale; dividing it by the acceleration to form a phase
    // duration multiplies that error by the reciprocal of the acceleration, and
    // the position the two phases have to agree on multiplies it by the velocity
    // once more. At a large velocity over a small acceleration this term is what
    // the two closed-form branches can differ by where they meet, and it is a
    // conditioning limit of the parametrization rather than a discontinuity in
    // the profile.
    const double seam_conditioning = std::numeric_limits<double>::epsilon() * rounding_op_margin
        * v_ceiling * v_ceiling / a_eff;

    // Sample-time floor: advancing through num_samples discrete sample times
    // carries its own fixed rounding floor of about one part in num_samples of
    // machine epsilon once divided back out by dt in a finite difference,
    // independent of which physical limit (v_max/a_eff) is being checked; it
    // only matters once that limit is itself very small.
    const double time_resolution_floor
        = std::numeric_limits<double>::epsilon() * rounding_op_margin * static_cast<double>(num_samples);

    double prev_q = 0.0;
    double prev_v = 0.0;
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
        // The velocity check adds an a_eff*dt physical term: when the
        // acceleration-to-cruise phase transition falls within a single
        // dt-sized sampling step (coarse sampling relative to that phase's
        // own, possibly tiny, duration), the returned v can differ from the
        // clean v_max/0 value by up to one step's worth of acceleration.
        const double v_env_tol = a_eff * dt
            + std::numeric_limits<double>::epsilon() * rounding_op_margin * v_max * static_cast<double>(num_samples)
            + time_resolution_floor;
        const double a_env_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin * a_eff * static_cast<double>(num_samples)
            + time_resolution_floor;
        if(std::abs(v) > v_max + v_env_tol)
            abort();
        if(std::abs(a) > a_eff + a_env_tol)
            abort();

        if(have_prev)
        {
            // Continuity: position cannot move faster than the kinematic bound
            // set by bounded velocity and acceleration over one sample step.
            const double continuity_bound = v_max * dt + 0.5 * a_eff * dt * dt;
            // The deceleration branch forms its position by subtracting the
            // distance still to run from the commanded displacement, so its
            // rounding lives at the scale of that displacement rather than at
            // the scale of the position value it returns; near the end of the
            // move those two differ by every digit.
            const double continuity_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin
                    * (std::abs(q) + std::abs(prev_q) + abs_h + continuity_bound)
                + seam_conditioning;
            if(std::abs(q - prev_q) > continuity_bound + continuity_tol)
                abort();

            // Finite-difference velocity from the position trace: within a
            // single accel/cruise/decel phase, position is quadratic in t so
            // the central-like forward difference matches the true derivative
            // up to a term bounded by one acceleration-phase jump (at most
            // a_eff) times the step; a small eps floor absorbs rounding. A
            // sample pair straddling a phase boundary is evaluated by two
            // distinct closed-form branches that meet exactly only in real
            // arithmetic, so an extra term scaled by the position magnitude in
            // play (the same rounding-op margin, divided back out by dt to
            // land in velocity units) absorbs that seam's rounding.
            const double v_fd = (q - prev_q) / dt;
            const double seam_v_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin
                    * (std::abs(q) + std::abs(prev_q) + abs_h) / dt
                + seam_conditioning / dt;
            const double v_fd_tol = a_eff * dt
                + std::numeric_limits<double>::epsilon() * rounding_op_margin * v_max * static_cast<double>(num_samples)
                + time_resolution_floor + seam_v_tol;
            if(std::abs(v_fd) > v_max + v_fd_tol)
                abort();

            // Finite-difference acceleration from the returned velocity trace:
            // velocity is piecewise linear in t (constant acceleration per
            // phase), so this difference is a convex combination of at most
            // two phases' accelerations and is bounded by a_eff up to
            // rounding, plus the same phase-boundary seam term one derivative
            // order up (scaled by the velocity magnitude in play).
            const double a_fd = (v - prev_v) / dt;
            const double seam_a_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin
                * (std::abs(v) + std::abs(prev_v)) / dt;
            const double a_fd_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin * a_eff * static_cast<double>(num_samples)
                + time_resolution_floor + seam_a_tol;
            if(std::abs(a_fd) > a_eff + a_fd_tol)
                abort();
        }

        prev_q = q;
        prev_v = v;
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

    ctrlpp::trapezoidal_trajectory<double> retimed = traj;
    const auto rescaled = retimed.rescale_to(T_new);
    if(!rescaled.has_value())
    {
        // A typed rejection is a correct outcome for a duration the command and
        // the limits cannot realize together, so it is not a finding. What would
        // be a finding is a rejection that left the profile changed anyway.
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
        * (static_cast<double>(area.panels) * area_scale + a_eff * T_scaled * area.back_time_exposure);
    if(!std::isfinite(displacement_tol))
        return 0;
    if(std::abs(area.integral - h_signed) > displacement_tol)
        abort();

    // Boundary velocities, a step inside each end. The step is the profile's own
    // first (respectively last) nonempty segment, and the excursion across it
    // cannot exceed the effective acceleration times its length.
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

    const double v_tol = std::numeric_limits<double>::epsilon() * rounding_op_margin
        * (v_ceiling + a_eff * T_scaled);
    if(!std::isfinite(v_tol))
        return 0;

    if(first_length > 0.0)
    {
        const double v_start = retimed.evaluate(first_length).velocity(0);
        if(!std::isfinite(v_start))
            return 0;
        if(std::abs(v_start - v0) > a_eff * first_length + v_tol)
            abort();
    }
    if(last_length > 0.0)
    {
        const double v_end = retimed.evaluate(T_scaled - last_length).velocity(0);
        if(!std::isfinite(v_end))
            return 0;
        if(std::abs(v_end - v1) > a_eff * last_length + v_tol)
            abort();
    }

    // Kinematic limits over the retimed profile, on a grid laid over the
    // requested duration rather than the realized one.
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
           > a_eff + std::numeric_limits<double>::epsilon() * rounding_op_margin * a_eff)
            abort();
    }

    // Realized duration. Thirty-two chained roundings stand behind it: twelve
    // form the cruise velocity (the linear coefficient, the constant, the
    // discriminant, its square root, and the cancellation-free root selection),
    // fifteen form the duration recomputed from it, and the remaining five cover
    // the domain clamps above and the multiplication that formed the request.
    // The solve's own conditioning multiplies all of them, because each shape
    // reaches its root through a difference of two nearly equal quantities.
    constexpr double duration_rounding_ops = 32.0;
    const double sigma = (h_signed >= 0.0) ? 1.0 : -1.0;
    const double conditioning = solve_conditioning(sigma * retimed.peak_velocity(), a_eff, abs_h,
                                                   sigma * v0, sigma * v1, T_new);
    const double duration_tol
        = std::numeric_limits<double>::epsilon() * duration_rounding_ops * T_new * conditioning;
    if(std::isfinite(duration_tol) && std::abs(T_scaled - T_new) > duration_tol)
        abort();

    return 0;
}
