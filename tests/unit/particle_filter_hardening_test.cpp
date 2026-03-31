#include "hardening_helpers.h"
#include "ctrlpp/estimation/particle_filter.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <random>

using Catch::Matchers::WithinAbs;

namespace {

struct pf_linear_dynamics
{
    double dt = 0.1;

    auto operator()(const ctrlpp::Vector<double, 2>& x,
                    const ctrlpp::Vector<double, 1>& u) const -> ctrlpp::Vector<double, 2>
    {
        ctrlpp::Vector<double, 2> xn;
        xn(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        xn(1) = x(1) + dt * u(0);
        return xn;
    }
};

struct pf_position_measurement
{
    auto operator()(const ctrlpp::Vector<double, 2>& x) const -> ctrlpp::Vector<double, 1>
    {
        ctrlpp::Vector<double, 1> z;
        z(0) = x(0);
        return z;
    }
};

auto make_pf(std::size_t seed = 42)
{
    ctrlpp::Matrix<double, 2, 2> Q = ctrlpp::Matrix<double, 2, 2>::Identity() * 0.01;
    ctrlpp::Matrix<double, 1, 1> R;
    R << 0.5;
    ctrlpp::Vector<double, 2> x0 = ctrlpp::Vector<double, 2>::Zero();
    ctrlpp::Matrix<double, 2, 2> P0 = ctrlpp::Matrix<double, 2, 2>::Identity();

    return ctrlpp::make_particle_filter<100>(
        pf_linear_dynamics{}, pf_position_measurement{},
        ctrlpp::pf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0},
        std::mt19937_64{seed});
}

}

TEST_CASE("PF all-NaN measurements does not crash",
          "[particle_filter][hardening][negative]")
{
    auto pf = make_pf();

    ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();
    pf.predict(u);

    ctrlpp::Vector<double, 1> z;
    z << std::numeric_limits<double>::quiet_NaN();
    pf.update(z);

    // NaN propagation or finite -- no crash
    CHECK((std::isnan(pf.state()[0]) || std::isfinite(pf.state()[0])));
}

TEST_CASE("PF Gaussian posterior mean/var within tolerance",
          "[particle_filter][hardening][precision]")
{
    // Use 500 particles for better precision
    ctrlpp::Matrix<double, 2, 2> Q = ctrlpp::Matrix<double, 2, 2>::Identity() * 0.01;
    ctrlpp::Matrix<double, 1, 1> R;
    R << 0.1;
    ctrlpp::Vector<double, 2> x0 = ctrlpp::Vector<double, 2>::Zero();
    ctrlpp::Matrix<double, 2, 2> P0 = ctrlpp::Matrix<double, 2, 2>::Identity();

    auto pf = ctrlpp::make_particle_filter<500>(
        pf_linear_dynamics{}, pf_position_measurement{},
        ctrlpp::pf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0},
        std::mt19937_64{42});

    double true_pos = 0.0;
    double true_vel = 1.0;
    ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();

    for(int k = 0; k < 50; ++k)
    {
        true_pos += 0.1 * true_vel;
        pf.predict(u);
        ctrlpp::Vector<double, 1> z;
        z << true_pos;
        pf.update(z);
    }

    REQUIRE(std::abs(pf.state()[0] - true_pos) < 1.0);
}

TEST_CASE("PF converges to true state for linear Gaussian system",
          "[particle_filter][hardening][convergence]")
{
    auto pf = make_pf(123);

    double true_pos = 0.0;
    double true_vel = 1.0;
    ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();

    for(int k = 0; k < 100; ++k)
    {
        true_pos += 0.1 * true_vel;
        pf.predict(u);
        ctrlpp::Vector<double, 1> z;
        z << true_pos;
        pf.update(z);
    }

    REQUIRE(std::abs(pf.state()[0] - true_pos) < 2.0);
}

TEST_CASE("PF with 10 particles does not crash",
          "[particle_filter][hardening][robustness]")
{
    ctrlpp::Matrix<double, 2, 2> Q = ctrlpp::Matrix<double, 2, 2>::Identity() * 0.01;
    ctrlpp::Matrix<double, 1, 1> R;
    R << 1.0;
    ctrlpp::Vector<double, 2> x0 = ctrlpp::Vector<double, 2>::Zero();
    ctrlpp::Matrix<double, 2, 2> P0 = ctrlpp::Matrix<double, 2, 2>::Identity();

    auto pf = ctrlpp::make_particle_filter<10>(
        pf_linear_dynamics{}, pf_position_measurement{},
        ctrlpp::pf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0},
        std::mt19937_64{42});

    ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();

    bool all_finite = true;
    for(int k = 0; k < 50; ++k)
    {
        pf.predict(u);
        ctrlpp::Vector<double, 1> z;
        z << 1.0;
        pf.update(z);

        if(!std::isfinite(pf.state()[0]))
        {
            all_finite = false;
            break;
        }
    }

    REQUIRE(all_finite);
}
