/// @file ctrlpp_ukf_01_pendulum.cpp
/// @brief Demonstrates UKF tracking a damped pendulum with noisy angle-only measurements.
///
/// Same pendulum system as the EKF example (state = [theta, omega], input = [torque],
/// measurement = [theta]) but using UKF -- no Jacobians needed. The dynamics and
/// measurement are plain callables satisfying dynamics_model and measurement_model.

#include "ctrlpp/estimation/ukf.h"

#include <cmath>
#include <cstddef>
#include <numbers>
#include <iostream>

namespace {

// Physical parameters
constexpr double g = 9.81;
constexpr double l = 1.0;
constexpr double b_damp = 0.1;
constexpr double m = 1.0;
constexpr double dt = 0.01;
constexpr std::size_t n_steps = 500;

// Process and measurement noise amplitudes
constexpr double process_noise_amp = 0.001;
constexpr double meas_noise_amp = 0.05;

/// Pendulum dynamics satisfying dynamics_model (but NOT differentiable_dynamics).
/// Forward-Euler discretisation: x_{k+1} = f(x_k, u_k).
/// No Jacobians -- UKF propagates sigma points through the nonlinear model directly.
struct pendulum_dynamics {
    auto operator()(const ctrlpp::Vector<double, 2>& x,
                    const ctrlpp::Vector<double, 1>& u) const -> ctrlpp::Vector<double, 2>
    {
        double theta = x(0);
        double omega = x(1);
        double tau = u(0);
        ctrlpp::Vector<double, 2> x_next;
        x_next(0) = theta + omega * dt;
        x_next(1) = omega + (-g / l * std::sin(theta) - b_damp * omega + tau / (m * l * l)) * dt;
        return x_next;
    }
};

/// Angle-only measurement satisfying measurement_model (no Jacobian).
struct angle_measurement {
    auto operator()(const ctrlpp::Vector<double, 2>& x) const -> ctrlpp::Vector<double, 1>
    {
        ctrlpp::Vector<double, 1> z;
        z(0) = x(0);
        return z;
    }
};

/// Simple deterministic pseudo-noise for reproducibility (linear congruential).
struct lcg_noise {
    std::uint32_t state;

    explicit lcg_noise(std::uint32_t seed) : state{seed} {}

    auto next() -> double
    {
        state = state * 1664525u + 1013904223u;
        // Map to [-1, 1]
        return static_cast<double>(static_cast<std::int32_t>(state)) / 2147483648.0;
    }
};

}

int main()
{
    pendulum_dynamics dyn;
    angle_measurement meas;

    ctrlpp::Matrix<double, 2, 2> Q = ctrlpp::Matrix<double, 2, 2>::Identity() * (process_noise_amp * process_noise_amp);
    ctrlpp::Matrix<double, 1, 1> R;
    R << meas_noise_amp * meas_noise_amp;
    ctrlpp::Vector<double, 2> x0 = ctrlpp::Vector<double, 2>::Zero();
    ctrlpp::Matrix<double, 2, 2> P0 = ctrlpp::Matrix<double, 2, 2>::Identity() * 1.0;

    ctrlpp::ukf filter(dyn, meas,
        ctrlpp::ukf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0});

    // True initial state: pendulum at 45 degrees, at rest
    ctrlpp::Vector<double, 2> x_true;
    x_true << std::numbers::pi / 4.0, 0.0;

    lcg_noise noise_gen{42};

    // CSV header
    std::cout << "step,true_theta,true_omega,est_theta,est_omega,P00,P11\n";

    for (std::size_t k = 0; k < n_steps; ++k) {
        // Zero torque input (free swing)
        ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();

        // Propagate true state with process noise
        double theta = x_true(0);
        double omega = x_true(1);
        x_true(0) = theta + omega * dt + process_noise_amp * noise_gen.next();
        x_true(1) = omega + (-g / l * std::sin(theta) - b_damp * omega) * dt
                     + process_noise_amp * noise_gen.next();

        // UKF predict
        filter.predict(u);

        // Noisy measurement of angle
        double meas_theta = x_true(0) + meas_noise_amp * noise_gen.next();
        ctrlpp::Vector<double, 1> z;
        z << meas_theta;

        // UKF update
        filter.update(z);

        auto est = filter.state();
        auto P = filter.covariance();

        std::cout << k << ','
                  << x_true(0) << ',' << x_true(1) << ','
                  << est(0) << ',' << est(1) << ','
                  << P(0, 0) << ',' << P(1, 1) << '\n';
    }

    // Final estimation error
    auto est_final = filter.state();
    double err_theta = std::abs(x_true(0) - est_final(0));
    double err_omega = std::abs(x_true(1) - est_final(1));
    std::cerr << "Final error: theta=" << err_theta << " omega=" << err_omega << '\n';

    return 0;
}
