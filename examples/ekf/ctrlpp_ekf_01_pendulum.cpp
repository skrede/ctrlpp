// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./ctrlpp_ekf_01_pendulum' using 1:2 with lines title 'true theta', '' using 1:3 with lines title 'estimated theta'"
// Redirect: ./ctrlpp_ekf_01_pendulum > output.csv

#include "ctrlpp/estimation/ekf.h"

#include <cmath>
#include <cstdint>
#include <numbers>
#include <iostream>

namespace
{

constexpr double g = 9.81;
constexpr double l = 1.0;
constexpr double b_damp = 0.1;
constexpr double m = 1.0;
constexpr double dt = 0.01;
constexpr std::int32_t n_steps = 500;

constexpr double process_noise_amp = 0.001;
constexpr double meas_noise_amp = 0.05;

// Forward-Euler discretization. Defining jacobian_x and jacobian_u satisfies differentiable_dynamics,
// so the EKF uses these analytical Jacobians instead of numerical ones.
struct pendulum_dynamics
{
    ctrlpp::Vector<double, 2> operator()(const ctrlpp::Vector<double, 2>& x, const ctrlpp::Vector<double, 1>& u) const
    {
        double theta = x(0);
        double omega = x(1);
        double tau = u(0);
        ctrlpp::Vector<double, 2> x_next;
        x_next(0) = theta + omega * dt;
        x_next(1) = omega + (-g / l * std::sin(theta) - b_damp * omega + tau / (m * l * l)) * dt;
        return x_next;
    }

    ctrlpp::Matrix<double, 2, 2> jacobian_x(const ctrlpp::Vector<double, 2>& x, const ctrlpp::Vector<double, 1>& /*u*/) const
    {
        ctrlpp::Matrix<double, 2, 2> F;
        F << 1.0, dt, -g / l * std::cos(x(0)) * dt, 1.0 - b_damp * dt;
        return F;
    }

    ctrlpp::Matrix<double, 2, 1> jacobian_u(const ctrlpp::Vector<double, 2>& /*x*/, const ctrlpp::Vector<double, 1>& /*u*/) const
    {
        ctrlpp::Matrix<double, 2, 1> G;
        G << 0.0, dt / (m * l * l);
        return G;
    }
};

struct angle_measurement
{
    ctrlpp::Vector<double, 1> operator()(const ctrlpp::Vector<double, 2>& x) const
    {
        ctrlpp::Vector<double, 1> z;
        z(0) = x(0);
        return z;
    }

    ctrlpp::Matrix<double, 1, 2> jacobian(const ctrlpp::Vector<double, 2>& /*x*/) const
    {
        ctrlpp::Matrix<double, 1, 2> H;
        H << 1.0, 0.0;
        return H;
    }
};

// Linear congruential generator with the Numerical Recipes constants, deterministic so every run prints
// the same output.
struct lcg_noise
{
    std::uint32_t state;

    explicit lcg_noise(std::uint32_t seed)
        : state{seed}
    {
    }

    double next()
    {
        state = state * 1664525u + 1013904223u;
        // Maps the state onto [-1, 1).
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

    auto filter_result = ctrlpp::ekf<double, 2, 1, 1, pendulum_dynamics, angle_measurement>::create(dyn, meas, ctrlpp::ekf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0});
    if(!filter_result.has_value())
    {
        std::cerr << "invalid EKF configuration\n";
        return 1;
    }
    auto& filter = *filter_result;

    ctrlpp::Vector<double, 2> x_true;
    x_true << std::numbers::pi / 4.0, 0.0;

    lcg_noise noise_gen{42};

    std::cout << "# time,true_theta,est_theta,true_omega,est_omega,meas_theta,P00,P11\n";

    for(std::int32_t k = 0; k < n_steps; ++k)
    {
        ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();

        double theta = x_true(0);
        double omega = x_true(1);
        x_true(0) = theta + omega * dt + process_noise_amp * noise_gen.next();
        x_true(1) = omega + (-g / l * std::sin(theta) - b_damp * omega) * dt + process_noise_amp * noise_gen.next();

        filter.predict(u);

        double meas_theta = x_true(0) + meas_noise_amp * noise_gen.next();
        ctrlpp::Vector<double, 1> z;
        z << meas_theta;

        // update refuses a non-finite sample rather than let it destroy the estimate, and reports whether
        // it ran.
        if(const auto stepped = filter.update(z); !stepped)
        {
            std::cerr << "EKF rejected the measurement at step " << k << "\n";
            return 1;
        }

        auto est = filter.state();
        auto P = filter.covariance();

        std::cout << static_cast<double>(k) * dt << ',' << x_true(0) << ',' << est(0) << ',' << x_true(1) << ',' << est(1) << ',' << meas_theta << ',' << P(0, 0) << ',' << P(1, 1)
                  << '\n';
    }

    return 0;
}
