#include "hardening_helpers.h"
#include "ctrlpp/estimation/mekf.h"
#include "ctrlpp/lie/so3.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <cmath>
#include <limits>
#include <numbers>

using Catch::Matchers::WithinAbs;

namespace {

struct gravity_measurement
{
    auto operator()(const Eigen::Quaternion<double>& q,
                    const ctrlpp::Vector<double, 3>& /*b*/) const -> ctrlpp::Vector<double, 3>
    {
        return q.toRotationMatrix().transpose().col(2);
    }
};

auto quat_angle(const Eigen::Quaterniond& q1, const Eigen::Quaterniond& q2) -> double
{
    return 2.0 * std::acos(std::min(1.0, std::abs(q1.dot(q2))));
}

auto make_mekf()
{
    ctrlpp::mekf_config<double, 3, 3> cfg;
    cfg.Q *= 1e-4;
    cfg.R *= 0.01;
    cfg.dt = 0.01;
    return ctrlpp::mekf(gravity_measurement{}, cfg);
}

}

TEST_CASE("MEKF zero rotation noise covariance", "[mekf][hardening][negative]")
{
    ctrlpp::mekf_config<double, 3, 3> cfg;
    cfg.Q = ctrlpp::Matrix<double, 6, 6>::Zero();
    cfg.R *= 0.01;
    cfg.dt = 0.01;

    auto filter = ctrlpp::mekf(gravity_measurement{}, cfg);

    ctrlpp::Vector<double, 3> gyro = ctrlpp::Vector<double, 3>::Zero();
    filter.predict(gyro);

    ctrlpp::Vector<double, 3> z;
    z << 0.0, 0.0, 1.0;
    filter.update(z);

    CHECK(std::isfinite(filter.state()[0]));
}

TEST_CASE("MEKF NaN quaternion measurement", "[mekf][hardening][negative]")
{
    auto filter = make_mekf();

    ctrlpp::Vector<double, 3> gyro = ctrlpp::Vector<double, 3>::Zero();
    filter.predict(gyro);

    ctrlpp::Vector<double, 3> z;
    z << std::numeric_limits<double>::quiet_NaN(), 0.0, 1.0;
    filter.update(z);

    // NaN propagation or finite -- no crash
    CHECK((std::isnan(filter.state()[0]) || std::isfinite(filter.state()[0])));
}

TEST_CASE("MEKF covariance stays PD over 1000 steps", "[mekf][hardening][stability]")
{
    auto filter = make_mekf();

    ctrlpp::Vector<double, 3> gyro = ctrlpp::Vector<double, 3>::Zero();
    bool all_pd = true;

    for(int k = 0; k < 1000; ++k)
    {
        filter.predict(gyro);
        ctrlpp::Vector<double, 3> z;
        z << 0.0, 0.0, 1.0;
        filter.update(z);

        Eigen::SelfAdjointEigenSolver<ctrlpp::Matrix<double, 6, 6>> eigsolver(filter.covariance());
        for(int i = 0; i < 6; ++i)
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

TEST_CASE("MEKF attitude converges to true orientation", "[mekf][hardening][convergence]")
{
    // Start tilted
    auto q_true = Eigen::Quaterniond::Identity();
    auto q_init = Eigen::Quaterniond(Eigen::AngleAxisd(0.3, Eigen::Vector3d::UnitX()));

    ctrlpp::mekf_config<double, 3, 3> cfg;
    cfg.Q *= 1e-4;
    cfg.R *= 0.01;
    cfg.dt = 0.01;
    cfg.q0 = q_init;

    auto filter = ctrlpp::mekf(gravity_measurement{}, cfg);

    ctrlpp::Vector<double, 3> gyro = ctrlpp::Vector<double, 3>::Zero();

    for(int k = 0; k < 500; ++k)
    {
        filter.predict(gyro);

        // Gravity in body frame for true attitude (identity = [0,0,1])
        ctrlpp::Vector<double, 3> z;
        z << 0.0, 0.0, 1.0;
        filter.update(z);
    }

    auto q_est = filter.attitude();
    double angle_error = quat_angle(q_est, q_true);
    REQUIRE(angle_error < 0.1);
}

TEST_CASE("MEKF 180-degree rotation (near singularity)", "[mekf][hardening][robustness]")
{
    // Start nearly 180 degrees from truth
    auto q_init = Eigen::Quaterniond(Eigen::AngleAxisd(
        std::numbers::pi - 0.01, Eigen::Vector3d::UnitZ()));

    ctrlpp::mekf_config<double, 3, 3> cfg;
    cfg.Q *= 1e-4;
    cfg.R *= 0.1;
    cfg.dt = 0.01;
    cfg.q0 = q_init;

    auto filter = ctrlpp::mekf(gravity_measurement{}, cfg);

    ctrlpp::Vector<double, 3> gyro = ctrlpp::Vector<double, 3>::Zero();
    bool all_finite = true;

    for(int k = 0; k < 200; ++k)
    {
        filter.predict(gyro);
        ctrlpp::Vector<double, 3> z;
        z << 0.0, 0.0, 1.0;
        filter.update(z);

        if(!std::isfinite(filter.state()[0]))
        {
            all_finite = false;
            break;
        }
    }

    REQUIRE(all_finite);
}
