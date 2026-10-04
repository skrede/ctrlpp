#ifndef HPP_GUARD_BENCHMARKS_PROBLEMS_QP_MPC_INITIAL_STATES_H
#define HPP_GUARD_BENCHMARKS_PROBLEMS_QP_MPC_INITIAL_STATES_H

#include <Eigen/Dense>

#include <cmath>
#include <cstdint>

namespace ctrlpp::bench::problems::qp
{

inline constexpr double state_phase_step = 0.7;
inline constexpr double closed_loop_phase_step = 0.8;

inline auto fixed_initial_state(std::int32_t nx, double amplitude) -> Eigen::VectorXd
{
    Eigen::VectorXd x0(nx);
    for(std::int32_t i = 0; i < nx; ++i)
        x0(i) = amplitude * std::cos(state_phase_step * i);
    return x0;
}

/// A prescribed sequence rather than a simulated trajectory, so every solver is handed the same problems.
inline auto closed_loop_initial_state(std::int32_t nx, double amplitude, std::int32_t step) -> Eigen::VectorXd
{
    Eigen::VectorXd x0(nx);
    for(std::int32_t i = 0; i < nx; ++i)
        x0(i) = amplitude * std::cos(closed_loop_phase_step * step + state_phase_step * i);
    return x0;
}

/// The churning benchmark's integer hash (a murmur3-style finalizer over a multiplicative step and
/// state mix), with a seed term that vanishes at seed 0 so seed 0 reproduces that benchmark's sequence.
inline auto churn_unit(std::int32_t seed, std::int32_t step, std::int32_t i) -> double
{
    std::uint32_t h = static_cast<std::uint32_t>(step) * 2654435761u + static_cast<std::uint32_t>(i) * 40503u
                    + 0x9e3779b9u + static_cast<std::uint32_t>(seed) * 0x85ebca6bu;
    h ^= h >> 15;
    h *= 2246822519u;
    h ^= h >> 13;
    h *= 3266489917u;
    h ^= h >> 16;
    return static_cast<double>(h) / 4294967296.0;
}

inline auto churning_initial_state(std::int32_t nx, double amplitude, std::int32_t seed, std::int32_t step)
    -> Eigen::VectorXd
{
    Eigen::VectorXd x0(nx);
    for(std::int32_t i = 0; i < nx; ++i)
        x0(i) = amplitude * (2.0 * churn_unit(seed, step, i) - 1.0);
    return x0;
}

}

#endif
