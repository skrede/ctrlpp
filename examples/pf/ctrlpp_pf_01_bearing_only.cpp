/// @file ctrlpp_pf_01_bearing_only.cpp
/// @brief Demonstrates particle filter tracking a target from bearing-only measurements.
///
/// System: state = [px, py, vx, vy] (2D position + velocity), input = zero,
/// measurement = bearing angle from a known sensor position. The bearing-only
/// measurement creates a banana-shaped posterior that KF/EKF/UKF cannot represent
/// well, making this the classic particle filter showcase.

#include "ctrlpp/estimation/particle_filter.h"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <numbers>
#include <random>

namespace
{

constexpr double dt = 0.1;
constexpr std::size_t n_steps = 200;
constexpr std::size_t n_particles = 1000;

// Sensor position (known)
constexpr double sensor_x = 0.0;
constexpr double sensor_y = 0.0;

/// Constant velocity dynamics: x_{k+1} = A * x_k (no control input used).
struct cv_dynamics
{
    auto operator()(const ctrlpp::Vector<double, 4>& x, const ctrlpp::Vector<double, 1>& /*u*/) const -> ctrlpp::Vector<double, 4>
    {
        ctrlpp::Vector<double, 4> x_next;
        x_next(0) = x(0) + x(2) * dt;
        x_next(1) = x(1) + x(3) * dt;
        x_next(2) = x(2);
        x_next(3) = x(3);
        return x_next;
    }
};

/// Bearing measurement: theta = atan2(py - sy, px - sx).
struct bearing_measurement
{
    double sx;
    double sy;

    auto operator()(const ctrlpp::Vector<double, 4>& x) const -> ctrlpp::Vector<double, 1>
    {
        ctrlpp::Vector<double, 1> z;
        z(0) = std::atan2(x(1) - sy, x(0) - sx);
        return z;
    }
};

/// Simple deterministic pseudo-noise for reproducibility (linear congruential).
struct lcg_noise
{
    std::uint32_t state;

    explicit lcg_noise(std::uint32_t seed) : state{seed} {}

    auto next() -> double
    {
        state = state * 1664525u + 1013904223u;
        return static_cast<double>(static_cast<std::int32_t>(state)) / 2147483648.0;
    }
};

} // namespace

int main()
{
    cv_dynamics dyn;
    bearing_measurement meas{.sx = sensor_x, .sy = sensor_y};

    // Process noise: small velocity perturbations
    ctrlpp::Matrix<double, 4, 4> Q = ctrlpp::Matrix<double, 4, 4>::Zero();
    Q(0, 0) = 0.01; // position noise
    Q(1, 1) = 0.01;
    Q(2, 2) = 0.001; // velocity noise
    Q(3, 3) = 0.001;

    // Bearing measurement noise: ~0.1 rad std dev
    ctrlpp::Matrix<double, 1, 1> R;
    R << 0.01; // 0.1^2

    // Initial state estimate: target starts at (10, 5) moving (1, 0.5) m/s
    ctrlpp::Vector<double, 4> x0;
    x0 << 10.0, 5.0, 1.0, 0.5;

    // Fairly uncertain initial covariance
    ctrlpp::Matrix<double, 4, 4> P0 = ctrlpp::Matrix<double, 4, 4>::Identity() * 4.0;

    ctrlpp::pf_config<double, 4, 1, 1> config{.Q = Q, .R = R, .x0 = x0, .P0 = P0};

    auto filter = ctrlpp::make_particle_filter<n_particles>(dyn, meas, config, std::mt19937_64{42});

    // True initial state (offset from estimate to test convergence)
    ctrlpp::Vector<double, 4> x_true;
    x_true << 12.0, 4.0, 0.8, 0.6;

    lcg_noise noise_gen{123};

    // CSV header
    std::cout << "step,true_px,true_py,est_px,est_py,true_vx,true_vy,est_vx,est_vy\n";

    for(std::size_t k = 0; k < n_steps; ++k)
    {
        // Zero input (constant velocity model)
        ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();

        // Propagate true state with small process noise
        double px = x_true(0) + x_true(2) * dt + 0.1 * noise_gen.next();
        double py = x_true(1) + x_true(3) * dt + 0.1 * noise_gen.next();
        double vx = x_true(2) + 0.01 * noise_gen.next();
        double vy = x_true(3) + 0.01 * noise_gen.next();
        x_true << px, py, vx, vy;

        // PF predict
        filter.predict(u);

        // Noisy bearing measurement
        double true_bearing = std::atan2(x_true(1) - sensor_y, x_true(0) - sensor_x);
        double meas_bearing = true_bearing + 0.1 * noise_gen.next();

        // Wrap to [-pi, pi]
        while(meas_bearing > std::numbers::pi)
            meas_bearing -= 2.0 * std::numbers::pi;
        while(meas_bearing < -std::numbers::pi)
            meas_bearing += 2.0 * std::numbers::pi;

        ctrlpp::Vector<double, 1> z;
        z << meas_bearing;

        // PF update
        filter.update(z);

        auto est = filter.state();

        std::cout << k << ',' << x_true(0) << ',' << x_true(1) << ',' << est(0) << ',' << est(1) << ',' << x_true(2) << ',' << x_true(3) << ',' << est(2) << ',' << est(3) << '\n';
    }

    // Final position error
    auto est_final = filter.state();
    double err_pos = std::sqrt((x_true(0) - est_final(0)) * (x_true(0) - est_final(0)) + (x_true(1) - est_final(1)) * (x_true(1) - est_final(1)));
    std::cerr << "Final position error: " << err_pos << " m\n";

    return 0;
}
