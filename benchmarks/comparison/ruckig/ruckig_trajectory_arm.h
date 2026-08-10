#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_RUCKIG_RUCKIG_TRAJECTORY_ARM_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_RUCKIG_RUCKIG_TRAJECTORY_ARM_H

#include "trajectory_scan.h"

#include <ruckig/ruckig.hpp>

#include <array>
#include <cstddef>

namespace ctrlpp::bench
{

inline ruckig::InputParameter<1> build_input(motion_command const& cmd)
{
    ruckig::InputParameter<1> input;
    input.current_position = {cmd.q0};
    input.current_velocity = {0.0};
    input.current_acceleration = {0.0};
    input.target_position = {cmd.q1};
    input.target_velocity = {0.0};
    input.target_acceleration = {0.0};
    input.max_velocity = {cmd.v_max};
    input.max_acceleration = {cmd.a_max};
    input.max_jerk = {cmd.j_max};
    return input;
}

// The single-call Ruckig::update path is deliberately absent from this arm. It
// recalculates the profile, then samples it through Trajectory::at_time, then
// reads a steady clock twice to fill OutputParameter::calculation_duration, and it
// copies and compares a whole InputParameter to drive its own state machine. A row
// built on it charges this side for two operations plus its own instrumentation
// while charging the other side for one.
class ruckig_trajectory_arm
{
public:
    explicit ruckig_trajectory_arm(motion_command const& cmd)
        : m_otg{cmd.control_period}
        , m_input{build_input(cmd)}
        , m_trajectory{}
        , m_position{}
        , m_velocity{}
        , m_acceleration{}
        , m_section{}
        , m_period{cmd.control_period}
        , m_horizon{}
        , m_time{}
    {
    }

    ruckig::Result synthesize() { return m_otg.calculate(m_input, m_trajectory); }

    kinematic_state advance()
    {
        kinematic_state const sampled = state(m_time);
        wrap_time(m_time, m_period, m_horizon);
        return sampled;
    }

    kinematic_state state(double t)
    {
        m_trajectory.at_time(t, m_position, m_velocity, m_acceleration, m_section);
        return {m_position[0], m_velocity[0], m_acceleration[0]};
    }

    double duration() const { return m_trajectory.get_duration(); }

    void set_horizon(double horizon) { m_horizon = horizon; }

private:
    ruckig::Ruckig<1>         m_otg;
    ruckig::InputParameter<1> m_input;
    ruckig::Trajectory<1>     m_trajectory;
    std::array<double, 1>     m_position;
    std::array<double, 1>     m_velocity;
    std::array<double, 1>     m_acceleration;
    std::size_t               m_section;
    double                    m_period;
    double                    m_horizon;
    double                    m_time;
};

}

#endif
