#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_RUCKIG_TRAJECTORY_SCAN_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_RUCKIG_TRAJECTORY_SCAN_H

#include <cmath>
#include <cstdio>
#include <limits>
#include <cstdint>
#include <algorithm>

namespace ctrlpp::bench
{

struct motion_command
{
    double q0;
    double q1;
    double v_max;
    double a_max;
    double j_max;
    double control_period;
};

struct kinematic_state
{
    double q;
    double v;
    double a;
};

struct profile_scan
{
    double peak_v;
    double peak_a;
    double peak_j;
    double symmetry;
    double end_position;
};

struct comparison_figures
{
    profile_scan library;
    profile_scan rival;
    double       step;
    double       horizon;
    double       duration_gap;
    double       sampled_gap;
};

inline void wrap_time(double& t, double period, double horizon)
{
    t += period;
    if(t >= horizon)
        t = 0.0;
}

// The coarsest step at which a straight line drawn between two samples is
// indistinguishable at double resolution from the velocity trace it spans. That
// trace's second derivative is bounded by the jerk limit, so linear interpolation
// errs by at most j_max*step^2/8; equating that to the resolution of the velocity
// limit fixes the step and leaves no sample count to invent.
inline double scan_step(motion_command const& cmd)
{
    return std::sqrt(8.0 * cmd.v_max * std::numeric_limits<double>::epsilon() / cmd.j_max);
}

// A rest-to-rest time-optimal profile under limits symmetric in sign has an
// acceleration odd about its own midpoint, hence a velocity even about it and a
// position point-symmetric through it. Neither library encodes that: both build
// and evaluate the accelerating and decelerating halves through separate
// expressions, so what the identity leaves is the arm's own arithmetic error.
inline double symmetry_residual(motion_command const& cmd, kinematic_state const& near,
                                kinematic_state const& far)
{
    double const span = std::abs(cmd.q1 - cmd.q0);
    return std::max({std::abs(near.q + far.q - (cmd.q0 + cmd.q1)) / span,
                     std::abs(near.v - far.v) / cmd.v_max,
                     std::abs(near.a + far.a) / cmd.a_max});
}

// A uniform scan cannot see a peak falling between two samples, and the
// finite-difference jerk subtracts two nearly equal accelerations, so each margin
// is taken net of the resolution it was observed at: a_max*step on the velocity
// peak, j_max*step on the acceleration peak, 2*a_max*eps/step on the jerk. A
// time-optimal profile saturates its jerk limit, so that margin sits at zero by
// construction and only the subtraction can certify it.
inline double worst_net_margin(profile_scan const& scan, motion_command const& cmd, double step)
{
    double const eps = std::numeric_limits<double>::epsilon();
    return std::max({(scan.peak_v - cmd.v_max - cmd.a_max * step) / cmd.v_max,
                     (scan.peak_a - cmd.a_max - cmd.j_max * step) / cmd.a_max,
                     (scan.peak_j - cmd.j_max - 2.0 * cmd.a_max * eps / step) / cmd.j_max});
}

// Scanned over the first half of the horizon with each sample paired against its
// mirror, so one pass covers the whole profile and hands the symmetry residual the
// two instants it is defined on.
template <typename Arm>
profile_scan scan_arm(Arm& arm, motion_command const& cmd, double step)
{
    double const horizon = arm.duration();
    int64_t const steps = static_cast<int64_t>(horizon / (2.0 * step));
    profile_scan scan{};
    double previous_a = arm.state(0.0).a;
    for(int64_t i = 1; i <= steps; ++i)
    {
        kinematic_state const near = arm.state(static_cast<double>(i) * step);
        kinematic_state const far = arm.state(horizon - static_cast<double>(i) * step);
        scan.peak_v = std::max({scan.peak_v, std::abs(near.v), std::abs(far.v)});
        scan.peak_a = std::max({scan.peak_a, std::abs(near.a), std::abs(far.a)});
        scan.peak_j = std::max(scan.peak_j, std::abs(near.a - previous_a) / step);
        scan.symmetry = std::max(scan.symmetry, symmetry_residual(cmd, near, far));
        previous_a = near.a;
    }
    scan.end_position = std::abs(arm.state(std::nextafter(horizon, 0.0)).q - cmd.q1);
    return scan;
}

template <typename ArmA, typename ArmB>
double sampled_deviation(ArmA& lhs, ArmB& rhs, motion_command const& cmd, double step, double horizon)
{
    double const span = std::abs(cmd.q1 - cmd.q0);
    int64_t const steps = static_cast<int64_t>(horizon / step);
    double worst = 0.0;
    for(int64_t i = 0; i <= steps; ++i)
    {
        double const t = static_cast<double>(i) * step;
        kinematic_state const a = lhs.state(t);
        kinematic_state const b = rhs.state(t);
        worst = std::max({worst, std::abs(a.q - b.q) / span, std::abs(a.v - b.v) / cmd.v_max,
                          std::abs(a.a - b.a) / cmd.a_max});
    }
    return worst;
}

template <typename ArmA, typename ArmB>
comparison_figures measure(ArmA& library, ArmB& rival, motion_command const& cmd)
{
    double const step = scan_step(cmd);
    profile_scan const library_scan = scan_arm(library, cmd, step);
    profile_scan const rival_scan = scan_arm(rival, cmd, step);
    double const horizon = std::min(library.duration(), rival.duration());
    return {library_scan,
            rival_scan,
            step,
            horizon,
            std::abs(library.duration() - rival.duration()),
            sampled_deviation(library, rival, cmd, step, horizon)};
}

inline void report_scan(char const* label, profile_scan const& scan, motion_command const& cmd,
                        double step)
{
    double const margin = worst_net_margin(scan, cmd, step);
    std::printf("%s peak |v|=%.17g |a|=%.17g |j|=%.17g, worst net margin=%.17g\n", label,
                scan.peak_v, scan.peak_a, scan.peak_j, margin);
    if(margin > 0.0)
        std::printf("LIMIT VIOLATION: %s leaves a shared kinematic limit by more than the scan "
                    "resolution can explain; its synthesis row is withdrawn\n",
                    label);
}

}

#endif
