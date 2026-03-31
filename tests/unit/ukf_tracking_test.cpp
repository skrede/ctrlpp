#include "ctrlpp/estimation/ukf.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>
#include <numbers>

using namespace ctrlpp;

// ---------------------------------------------------------------------------
// Test dynamics: linear constant-velocity (double integrator)
// ---------------------------------------------------------------------------
struct ukf_linear_dynamics
{
    double dt = 0.1;

    auto operator()(const Vector<double, 2>& x, const Vector<double, 1>& u) const -> Vector<double, 2>
    {
        Vector<double, 2> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    }
};

struct ukf_position_measurement
{
    auto operator()(const Vector<double, 2>& x) const -> Vector<double, 1>
    {
        Vector<double, 1> z;
        z(0) = x(0);
        return z;
    }
};

// ---------------------------------------------------------------------------
// Test dynamics/measurement: pendulum (nonlinear)
// ---------------------------------------------------------------------------
struct ukf_pendulum_dynamics
{
    static constexpr double g = 9.81;
    static constexpr double l = 1.0;
    static constexpr double b = 0.1;
    static constexpr double m = 1.0;
    static constexpr double dt = 0.01;

    auto operator()(const Vector<double, 2>& x, const Vector<double, 1>& u) const -> Vector<double, 2>
    {
        double theta = x(0);
        double omega = x(1);
        double tau = u(0);
        Vector<double, 2> x_next;
        x_next(0) = theta + omega * dt;
        x_next(1) = omega + (-g / l * std::sin(theta) - b * omega + tau / (m * l * l)) * dt;
        return x_next;
    }
};

struct ukf_angle_measurement
{
    auto operator()(const Vector<double, 2>& x) const -> Vector<double, 1>
    {
        Vector<double, 1> z;
        z(0) = x(0);
        return z;
    }
};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("ukf tracks linear system")
{
    ukf_linear_dynamics dyn;
    ukf_position_measurement meas;

    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.01;
    Matrix<double, 1, 1> R;
    R << 1.0;
    Vector<double, 2> x0 = Vector<double, 2>::Zero();
    Matrix<double, 2, 2> P0 = Matrix<double, 2, 2>::Identity() * 10.0;

    ukf filter(dyn, meas, ukf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0});

    double true_pos = 0.0;
    double true_vel = 1.0;
    constexpr double dt = 0.1;

    double initial_trace = filter.covariance().trace();

    for(int i = 0; i < 50; ++i)
    {
        true_pos += true_vel * dt;

        Vector<double, 1> u = Vector<double, 1>::Zero();
        filter.predict(u);

        Vector<double, 1> z;
        z << true_pos + 0.1 * std::sin(static_cast<double>(i));
        filter.update(z);
    }

    auto est = filter.state();
    CHECK_THAT(est(0), Catch::Matchers::WithinAbs(true_pos, 0.5));
    CHECK_THAT(est(1), Catch::Matchers::WithinAbs(true_vel, 0.5));

    // Covariance should have converged (trace decreased)
    CHECK(filter.covariance().trace() < initial_trace);
}

TEST_CASE("ukf tracks nonlinear system")
{
    ukf_pendulum_dynamics dyn;
    ukf_angle_measurement meas;

    constexpr double dt = ukf_pendulum_dynamics::dt;
    constexpr double g = ukf_pendulum_dynamics::g;
    constexpr double l = ukf_pendulum_dynamics::l;
    constexpr double b = ukf_pendulum_dynamics::b;

    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.001;
    Matrix<double, 1, 1> R;
    R << 0.01;
    Vector<double, 2> x0;
    x0 << 0.1, 0.0;
    Matrix<double, 2, 2> P0 = Matrix<double, 2, 2>::Identity() * 1.0;

    ukf filter(dyn, meas, ukf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0});

    Vector<double, 2> x_true;
    x_true << std::numbers::pi / 4.0, 0.0;

    for(int i = 0; i < 200; ++i)
    {
        Vector<double, 1> u = Vector<double, 1>::Zero();

        Vector<double, 2> x_true_next;
        x_true_next(0) = x_true(0) + x_true(1) * dt;
        x_true_next(1) = x_true(1) + (-g / l * std::sin(x_true(0)) - b * x_true(1)) * dt;
        x_true = x_true_next;

        filter.predict(u);

        Vector<double, 1> z;
        z << x_true(0) + 0.05 * std::sin(static_cast<double>(i) * 0.7);
        filter.update(z);
    }

    auto est = filter.state();
    CHECK_THAT(est(0), Catch::Matchers::WithinAbs(x_true(0), 0.3));
    CHECK_THAT(est(1), Catch::Matchers::WithinAbs(x_true(1), 0.5));
}

TEST_CASE("ukf covariance remains PSD")
{
    ukf_linear_dynamics dyn;
    ukf_position_measurement meas;

    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.01;
    Matrix<double, 1, 1> R;
    R << 1.0;
    Vector<double, 2> x0 = Vector<double, 2>::Zero();
    Matrix<double, 2, 2> P0 = Matrix<double, 2, 2>::Identity() * 10.0;

    ukf filter(dyn, meas, ukf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0});

    for(int i = 0; i < 150; ++i)
    {
        Vector<double, 1> u = Vector<double, 1>::Zero();
        filter.predict(u);

        Vector<double, 1> z;
        z << static_cast<double>(i) * 0.1;
        filter.update(z);

        auto P = filter.covariance();
        CHECK((P - P.transpose()).norm() < 1e-10);
        Eigen::SelfAdjointEigenSolver<Matrix<double, 2, 2>> eigsolver(P, Eigen::EigenvaluesOnly);
        for(int j = 0; j < 2; ++j)
            CHECK(eigsolver.eigenvalues()(j) >= -1e-10);
    }
}
