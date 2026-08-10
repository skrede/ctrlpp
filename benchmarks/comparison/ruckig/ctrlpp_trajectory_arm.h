#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_RUCKIG_CTRLPP_TRAJECTORY_ARM_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_RUCKIG_CTRLPP_TRAJECTORY_ARM_H

#include "trajectory_scan.h"

#include "bench_construct.h"

#include "ctrlpp/expected.h"

#include "ctrlpp/trajectory/online_planner_3rd.h"

namespace ctrlpp::bench
{

class ctrlpp_trajectory_arm
{
public:
    explicit ctrlpp_trajectory_arm(motion_command const& cmd)
        : m_planner{built_or_exit(ctrlpp::online_planner_3rd<double>::create(
                                      {.v_max = cmd.v_max, .a_max = cmd.a_max, .j_max = cmd.j_max}),
                                  "ctrlpp::online_planner_3rd")}
        , m_command{cmd}
        , m_horizon{}
        , m_time{}
    {
    }

    ctrlpp::expected<void, trajectory_error> synthesize() { return m_planner.update(m_command.q1); }

    kinematic_state advance()
    {
        kinematic_state const sampled = state(m_time);
        wrap_time(m_time, m_command.control_period, m_horizon);
        return sampled;
    }

    kinematic_state state(double t)
    {
        auto const point = m_planner.sample(t);
        return {point.position[0], point.velocity[0], point.acceleration[0]};
    }

    double duration() const { return m_planner.diagnostics().planned_duration; }

    void set_horizon(double horizon) { m_horizon = horizon; }

    // update() replans from wherever sample() last left the planner, so a scan of
    // the profile has to be undone before the timed synthesis rows: without this,
    // every one of them would replan an axis already standing on its target.
    void rewind()
    {
        m_planner.reset(m_command.q0);
        m_time = 0.0;
    }

private:
    ctrlpp::online_planner_3rd<double> m_planner;
    motion_command                     m_command;
    double                             m_horizon;
    double                             m_time;
};

}

#endif
