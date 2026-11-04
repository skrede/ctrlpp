// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./ctrlpp_pf_01_bearing_only' using 1:2 with lines title 'true px', '' using 1:3 with lines title 'estimated px'"
// Redirect: ./ctrlpp_pf_01_bearing_only > output.csv
// A bearing-only measurement gives a banana-shaped posterior that a Kalman-type filter cannot represent,
// which is why this example uses a particle filter.

#include "ctrlpp/estimation/particle_filter.h"

#include <cmath>
#include <random>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <iostream>

namespace
{

constexpr double dt = 0.1;
constexpr std::int32_t n_steps = 200;
constexpr std::size_t n_particles = 1000;

constexpr double sensor_x = 0.0;
constexpr double sensor_y = 0.0;

// Constant-velocity model; the input is unused.
struct cv_dynamics
{
    ctrlpp::Vector<double, 4> operator()(const ctrlpp::Vector<double, 4>& x, const ctrlpp::Vector<double, 1>& /*u*/) const
    {
        ctrlpp::Vector<double, 4> x_next;
        x_next(0) = x(0) + x(2) * dt;
        x_next(1) = x(1) + x(3) * dt;
        x_next(2) = x(2);
        x_next(3) = x(3);
        return x_next;
    }
};

struct bearing_measurement
{
    double sx;
    double sy;

    ctrlpp::Vector<double, 1> operator()(const ctrlpp::Vector<double, 4>& x) const
    {
        ctrlpp::Vector<double, 1> z;
        z(0) = std::atan2(x(1) - sy, x(0) - sx);
        return z;
    }
};

// A bearing innovation formed as a plain difference is wrong across the +/-pi cut, where the default
// Gaussian likelihood treats adjacent angles as maximally wrong and collapses the filter whenever the target
// crosses it. This policy wraps the innovation with std::remainder before forming the Gaussian.
struct wrapped_bearing_likelihood
{
    double operator()(const ctrlpp::Vector<double, 1>& z, const ctrlpp::Vector<double, 1>& z_pred, const ctrlpp::Matrix<double, 1, 1>& R_inv, double log_det_2piR) const
    {
        double innov = std::remainder(z(0) - z_pred(0), 2.0 * std::numbers::pi);
        double mahal = innov * R_inv(0, 0) * innov;
        return -0.5 * mahal - 0.5 * log_det_2piR;
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
        return static_cast<double>(static_cast<std::int32_t>(state)) / 2147483648.0;
    }
};

const char* refusal_cause(ctrlpp::pf_update_error error)
{
    switch(error)
    {
        case ctrlpp::pf_update_error::non_finite_state:
            return "the particles are non-finite";
        case ctrlpp::pf_update_error::non_finite_measurement:
            return "the measurement is non-finite";
        case ctrlpp::pf_update_error::non_finite_result:
            return "the likelihood is non-finite";
    }
    return "unknown";
}

const char* health_name(ctrlpp::pf_health health)
{
    switch(health)
    {
        case ctrlpp::pf_health::ok:
            return "ok";
        case ctrlpp::pf_health::weights_reset:
            return "weights reset by a measurement no particle explained";
        case ctrlpp::pf_health::resampling_rolled_back:
            return "resampling rolled back after a non-finite roughening";
        case ctrlpp::pf_health::non_finite_estimate:
            return "non-finite estimate";
    }
    return "unknown";
}

}

int main()
{
    cv_dynamics dyn;
    bearing_measurement meas{.sx = sensor_x, .sy = sensor_y};

    ctrlpp::Matrix<double, 4, 4> Q = ctrlpp::Matrix<double, 4, 4>::Zero();
    Q(0, 0) = 0.01; // position noise
    Q(1, 1) = 0.01;
    Q(2, 2) = 0.001; // velocity noise
    Q(3, 3) = 0.001;

    ctrlpp::Matrix<double, 1, 1> R;
    R << 0.01; // (0.1 rad)^2

    ctrlpp::Vector<double, 4> x0;
    x0 << 10.0, 5.0, 1.0, 0.5;

    ctrlpp::Matrix<double, 4, 4> P0 = ctrlpp::Matrix<double, 4, 4>::Identity() * 4.0;

    ctrlpp::pf_config<double, 4, 1, 1> config{.Q = Q, .R = R, .x0 = x0, .P0 = P0};

    auto filter_result = ctrlpp::make_particle_filter<n_particles>(dyn, meas, config, std::mt19937_64{42}, wrapped_bearing_likelihood{});
    if(!filter_result.has_value())
    {
        std::cerr << "invalid particle filter configuration\n";
        return 1;
    }
    auto& filter = *filter_result;

    // Offset from the initial estimate, so the output shows the filter converging.
    ctrlpp::Vector<double, 4> x_true;
    x_true << 12.0, 4.0, 0.8, 0.6;

    lcg_noise noise_gen{123};

    std::cout << "# time,true_px,est_px,true_py,est_py,true_vx,est_vx,true_vy,est_vy\n";

    for(std::int32_t k = 0; k < n_steps; ++k)
    {
        ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();

        double px = x_true(0) + x_true(2) * dt + 0.1 * noise_gen.next();
        double py = x_true(1) + x_true(3) * dt + 0.1 * noise_gen.next();
        double vx = x_true(2) + 0.01 * noise_gen.next();
        double vy = x_true(3) + 0.01 * noise_gen.next();
        x_true << px, py, vx, vy;

        filter.predict(u);

        double true_bearing = std::atan2(x_true(1) - sensor_y, x_true(0) - sensor_x);
        double meas_bearing = true_bearing + 0.1 * noise_gen.next();

        while(meas_bearing > std::numbers::pi)
            meas_bearing -= 2.0 * std::numbers::pi;
        while(meas_bearing < -std::numbers::pi)
            meas_bearing += 2.0 * std::numbers::pi;

        ctrlpp::Vector<double, 1> z;
        z << meas_bearing;

        // A refused measurement leaves the particles, weights and generator untouched, so the filter carries
        // on with the next sample.
        if(const auto stepped = filter.update(z); !stepped)
        {
            std::cerr << "Measurement refused at step " << k << ": " << refusal_cause(stepped.error()) << '\n';
            continue;
        }

        auto est = filter.state();

        std::cout << static_cast<double>(k) * dt << ',' << x_true(0) << ',' << est(0) << ',' << x_true(1) << ',' << est(1) << ',' << x_true(2) << ',' << est(2) << ',' << x_true(3)
                  << ',' << est(3) << '\n';
    }

    // The final error is reported beside the spread the filter claims, since it is only interpretable
    // against it: the same error is healthy for a meter of spread and a failure for a centimeter.
    auto est_final = filter.state();
    double err_pos = std::sqrt((x_true(0) - est_final(0)) * (x_true(0) - est_final(0)) + (x_true(1) - est_final(1)) * (x_true(1) - est_final(1)));
    auto P_final = filter.covariance();
    std::cerr << "Final position error: " << err_pos << " m\n";
    std::cerr << "Reported position standard deviation: " << std::sqrt(P_final(0, 0)) << ", " << std::sqrt(P_final(1, 1)) << " m\n";
    std::cerr << "Filter health: " << health_name(filter.health()) << '\n';

    return 0;
}
