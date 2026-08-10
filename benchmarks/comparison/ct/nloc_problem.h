#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_CT_NLOC_PROBLEM_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_CT_NLOC_PROBLEM_H

#include <Eigen/Dense>

#include <cstddef>
#include <cstdint>

namespace ctrlpp::bench
{

constexpr std::size_t oscillator_state_dim = 2;
constexpr std::size_t oscillator_input_dim = 1;
constexpr int32_t oscillator_horizon = 20;
constexpr int32_t oscillator_sim_steps = 20;
constexpr double oscillator_sample_period = 0.05;
constexpr double oscillator_damping = 0.1;

using oscillator_state = Eigen::Matrix<double, int(oscillator_state_dim), 1>;
using oscillator_input = Eigen::Matrix<double, int(oscillator_input_dim), 1>;
using oscillator_state_weight = Eigen::Matrix<double, int(oscillator_state_dim), int(oscillator_state_dim)>;
using oscillator_input_weight = Eigen::Matrix<double, int(oscillator_input_dim), int(oscillator_input_dim)>;

inline oscillator_state_weight oscillator_dynamics_matrix()
{
    oscillator_state_weight A = oscillator_state_weight::Zero();
    A(0, 1) = 1.0;
    A(1, 0) = -1.0;
    A(1, 1) = -oscillator_damping;
    return A;
}

inline Eigen::Matrix<double, int(oscillator_state_dim), int(oscillator_input_dim)> oscillator_input_matrix()
{
    Eigen::Matrix<double, int(oscillator_state_dim), int(oscillator_input_dim)> B =
        Eigen::Matrix<double, int(oscillator_state_dim), int(oscillator_input_dim)>::Zero();
    B(1, 0) = 1.0;
    return B;
}

// The explicit Euler pair the competitor's forward-Euler discretization forms
// from the continuous system, so both arms advance the identical plant.
inline oscillator_state_weight oscillator_transition_matrix()
{
    return oscillator_state_weight::Identity() + oscillator_sample_period * oscillator_dynamics_matrix();
}

inline Eigen::Matrix<double, int(oscillator_state_dim), int(oscillator_input_dim)> oscillator_input_gain()
{
    return oscillator_sample_period * oscillator_input_matrix();
}

inline oscillator_state oscillator_step(const oscillator_state& x, const oscillator_input& u)
{
    return oscillator_transition_matrix() * x + oscillator_input_gain() * u;
}

inline oscillator_state oscillator_initial_state()
{
    return oscillator_state{1.0, 0.0};
}

struct loop_outcome
{
    double cost;
    double first_input;
};

// The one objective both arms are scored by. It is the sampled integral of the
// running cost, so it is the competitor's own continuous cost weighting and the
// weights handed to this repository's controller are scaled to match it.
inline double oscillator_stage_cost(const oscillator_state& x, const oscillator_input& u)
{
    return oscillator_sample_period * (x.squaredNorm() + u.squaredNorm());
}

}

#endif
