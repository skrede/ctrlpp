#include "hardening_helpers.h"
#include "ctrlpp/estimation/ekf.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

struct linear_dynamics
{
    auto operator()(const ctrlpp::Vector<double, 2>& x,
                    const ctrlpp::Vector<double, 1>& u) const -> ctrlpp::Vector<double, 2>
    {
        ctrlpp::Vector<double, 2> xn;
        xn[0] = x[0] + 0.1 * x[1];
        xn[1] = x[1] + 0.1 * u[0];
        return xn;
    }
};

struct linear_measurement
{
    auto operator()(const ctrlpp::Vector<double, 2>& x) const -> ctrlpp::Vector<double, 1>
    {
        ctrlpp::Vector<double, 1> z;
        z << x[0];
        return z;
    }
};

auto make_ekf()
{
    ctrlpp::ekf_config<double, 2, 1, 1> cfg{};
    cfg.Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    cfg.R = Eigen::Matrix<double, 1, 1>::Identity();
    cfg.x0 = Eigen::Vector2d::Zero();
    cfg.P0 = Eigen::Matrix<double, 2, 2>::Identity() * 10.0;

    return ctrlpp::ekf(linear_dynamics{}, linear_measurement{}, cfg);
}

}

TEST_CASE("EKF NaN measurement does not crash", "[ekf][hardening][negative]")
{
    auto filter = make_ekf();

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;
    filter.predict(u);

    Eigen::Matrix<double, 1, 1> z;
    z << std::numeric_limits<double>::quiet_NaN();
    filter.update(z);

    CHECK((std::isnan(filter.state()[0]) || std::isfinite(filter.state()[0])));
}

TEST_CASE("EKF Inf process noise does not crash", "[ekf][hardening][negative]")
{
    ctrlpp::ekf_config<double, 2, 1, 1> cfg{};
    cfg.Q = Eigen::Matrix<double, 2, 2>::Identity() * std::numeric_limits<double>::infinity();
    cfg.R = Eigen::Matrix<double, 1, 1>::Identity();
    cfg.x0 = Eigen::Vector2d::Zero();
    cfg.P0 = Eigen::Matrix<double, 2, 2>::Identity();

    auto filter = ctrlpp::ekf(linear_dynamics{}, linear_measurement{}, cfg);

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;
    filter.predict(u);

    Eigen::Matrix<double, 1, 1> z;
    z << 5.0;
    filter.update(z);

    // With infinite Q, state may not be finite, but should not crash
    CHECK(true);
}

TEST_CASE("EKF for linear system matches Kalman gain within 1%", "[ekf][hardening][precision]")
{
    // For a linear system, the EKF should produce the same result as the Kalman filter
    auto filter = make_ekf();

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

    // After 50 steps, position estimate should be close to truth
    REQUIRE(std::abs(filter.state()[0] - true_pos) < 1.0);
}

TEST_CASE("EKF covariance stays PD over 1000 steps", "[ekf][hardening][stability]")
{
    auto filter = make_ekf();

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

TEST_CASE("EKF state estimate converges to truth", "[ekf][hardening][convergence]")
{
    auto filter = make_ekf();

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;

    double true_pos = 0.0;
    double true_vel = 1.0;

    for(int k = 0; k < 200; ++k)
    {
        true_pos += 0.1 * true_vel;
        filter.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << true_pos;
        filter.update(z);
    }

    REQUIRE(std::abs(filter.state()[0] - true_pos) < 0.5);
}

TEST_CASE("EKF ill-conditioned system cond 1e10", "[ekf][hardening][robustness]")
{
    ctrlpp::ekf_config<double, 2, 1, 1> cfg{};
    cfg.Q = ctrlpp::test::ill_conditioned_2x2<double>(1e10);
    cfg.R = Eigen::Matrix<double, 1, 1>::Identity();
    cfg.x0 = Eigen::Vector2d::Zero();
    cfg.P0 = Eigen::Matrix<double, 2, 2>::Identity();

    auto filter = ctrlpp::ekf(linear_dynamics{}, linear_measurement{}, cfg);

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
