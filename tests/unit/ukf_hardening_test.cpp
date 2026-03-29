#include "hardening_helpers.h"
#include "ctrlpp/estimation/ukf.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;

namespace {

struct ukf_linear_dynamics
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

struct ukf_position_measurement
{
    auto operator()(const ctrlpp::Vector<double, 2>& x) const -> ctrlpp::Vector<double, 1>
    {
        ctrlpp::Vector<double, 1> z;
        z(0) = x(0);
        return z;
    }
};

auto make_ukf()
{
    ctrlpp::ukf_config<double, 2, 1, 1> cfg{};
    cfg.Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    cfg.R = Eigen::Matrix<double, 1, 1>::Identity();
    cfg.x0 = Eigen::Vector2d::Zero();
    cfg.P0 = Eigen::Matrix<double, 2, 2>::Identity() * 10.0;
    return ctrlpp::ukf(ukf_linear_dynamics{}, ukf_position_measurement{}, cfg);
}

}

TEST_CASE("UKF NaN measurement does not crash", "[ukf][hardening][negative]")
{
    auto filter = make_ukf();

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;
    filter.predict(u);

    Eigen::Matrix<double, 1, 1> z;
    z << std::numeric_limits<double>::quiet_NaN();
    filter.update(z);

    CHECK((std::isnan(filter.state()[0]) || std::isfinite(filter.state()[0])));
}

TEST_CASE("UKF for linear system matches Kalman output", "[ukf][hardening][precision]")
{
    auto filter = make_ukf();

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;

    double true_pos = 0.0;
    double true_vel = 1.0;

    for(int k = 0; k < 50; ++k)
    {
        true_pos += 0.1 * true_vel;
        filter.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << true_pos;
        filter.update(z);
    }

    REQUIRE(std::abs(filter.state()[0] - true_pos) < 1.0);
}

TEST_CASE("UKF covariance stays PD over 1000 steps", "[ukf][hardening][stability]")
{
    auto filter = make_ukf();

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;
    bool all_pd = true;

    for(int k = 0; k < 1000; ++k)
    {
        filter.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << 5.0 + 0.01 * k;
        filter.update(z);

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> eigsolver(filter.covariance());
        for(int i = 0; i < 2; ++i)
        {
            if(eigsolver.eigenvalues()(i) < -1e-10)
            {
                all_pd = false;
                break;
            }
        }
        if(!all_pd)
            break;
    }

    REQUIRE(all_pd);
}

TEST_CASE("UKF tracks nonlinear system (quadratic dynamics)", "[ukf][hardening][convergence]")
{
    struct quadratic_dynamics
    {
        auto operator()(const ctrlpp::Vector<double, 2>& x,
                        const ctrlpp::Vector<double, 1>& /*u*/) const -> ctrlpp::Vector<double, 2>
        {
            ctrlpp::Vector<double, 2> xn;
            xn(0) = 0.95 * x(0) + 0.05 * x(1);
            xn(1) = 0.95 * x(1);
            return xn;
        }
    };

    ctrlpp::ukf_config<double, 2, 1, 1> cfg{};
    cfg.Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    cfg.R = Eigen::Matrix<double, 1, 1>::Identity() * 0.1;
    cfg.x0 = Eigen::Vector2d::Zero();
    cfg.P0 = Eigen::Matrix<double, 2, 2>::Identity() * 10.0;

    auto filter = ctrlpp::ukf(quadratic_dynamics{}, ukf_position_measurement{}, cfg);

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;

    // True state converges to zero -- filter should track
    ctrlpp::Vector<double, 2> x_true;
    x_true << 5.0, 1.0;

    for(int k = 0; k < 200; ++k)
    {
        x_true(0) = 0.95 * x_true(0) + 0.05 * x_true(1);
        x_true(1) = 0.95 * x_true(1);

        filter.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << x_true(0);
        filter.update(z);
    }

    REQUIRE(std::abs(filter.state()[0] - x_true(0)) < 0.5);
}

TEST_CASE("UKF near-singular P0", "[ukf][hardening][robustness]")
{
    ctrlpp::ukf_config<double, 2, 1, 1> cfg{};
    cfg.Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    cfg.R = Eigen::Matrix<double, 1, 1>::Identity();
    cfg.x0 = Eigen::Vector2d::Zero();
    cfg.P0 = ctrlpp::test::ill_conditioned_2x2<double>(1e10);

    auto filter = ctrlpp::ukf(ukf_linear_dynamics{}, ukf_position_measurement{}, cfg);

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;

    bool all_finite = true;
    for(int k = 0; k < 100; ++k)
    {
        filter.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << 1.0;
        filter.update(z);

        if(!std::isfinite(filter.state()[0]) || !std::isfinite(filter.state()[1]))
        {
            all_finite = false;
            break;
        }
    }

    REQUIRE(all_finite);
}
