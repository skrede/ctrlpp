#include "ctrlpp/estimation/ekf.h"
#include "ctrlpp/estimation/observer_policy.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>

using namespace ctrlpp;

// ---------------------------------------------------------------------------
// Test dynamics: linear constant-velocity model with analytical Jacobians
// ---------------------------------------------------------------------------
struct linear_dynamics
{
    double dt = 0.1;

    auto operator()(const Vector<double, 2>& x, const Vector<double, 1>& u) const -> Vector<double, 2>
    {
        Vector<double, 2> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    }

    auto jacobian_x(const Vector<double, 2>& /*x*/, const Vector<double, 1>& /*u*/) const -> Matrix<double, 2, 2>
    {
        Matrix<double, 2, 2> F;
        F << 1.0, dt, 0.0, 1.0;
        return F;
    }

    auto jacobian_u(const Vector<double, 2>& /*x*/, const Vector<double, 1>& /*u*/) const -> Matrix<double, 2, 1>
    {
        Matrix<double, 2, 1> G;
        G << 0.5 * dt * dt, dt;
        return G;
    }
};

// ---------------------------------------------------------------------------
// Test measurement: position observation with analytical Jacobian
// ---------------------------------------------------------------------------
struct position_measurement
{
    auto operator()(const Vector<double, 2>& x) const -> Vector<double, 1>
    {
        Vector<double, 1> z;
        z(0) = x(0);
        return z;
    }

    auto jacobian(const Vector<double, 2>& /*x*/) const -> Matrix<double, 1, 2>
    {
        Matrix<double, 1, 2> H;
        H << 1.0, 0.0;
        return H;
    }
};

// ---------------------------------------------------------------------------
// Concept satisfaction: lambda-based EKF (numerical Jacobians)
// ---------------------------------------------------------------------------
using lambda_dyn_t = decltype([](const Vector<double, 2>& x,
                                 const Vector<double, 1>& u) -> Vector<double, 2> {
    constexpr double dt = 0.1;
    Vector<double, 2> x_next;
    x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
    x_next(1) = x(1) + dt * u(0);
    return x_next;
});

using lambda_meas_t = decltype([](const Vector<double, 2>& x) -> Vector<double, 1> {
    Vector<double, 1> z;
    z(0) = x(0);
    return z;
});

static_assert(ObserverPolicy<ekf<double, 2, 1, 1, lambda_dyn_t, lambda_meas_t>>);
static_assert(CovarianceObserver<ekf<double, 2, 1, 1, lambda_dyn_t, lambda_meas_t>>);

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("ekf covariance stays symmetric and PSD over 100+ cycles")
{
    linear_dynamics dyn;
    position_measurement meas;

    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.01;
    Matrix<double, 1, 1> R;
    R << 1.0;
    Vector<double, 2> x0 = Vector<double, 2>::Zero();
    Matrix<double, 2, 2> P0 = Matrix<double, 2, 2>::Identity() * 10.0;

    ekf filter(dyn, meas, ekf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0});

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

TEST_CASE("ekf NEES is finite and positive")
{
    linear_dynamics dyn;
    position_measurement meas;

    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.01;
    Matrix<double, 1, 1> R;
    R << 1.0;
    Vector<double, 2> x0 = Vector<double, 2>::Zero();
    Matrix<double, 2, 2> P0 = Matrix<double, 2, 2>::Identity() * 10.0;

    ekf filter(dyn, meas, ekf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0});

    Vector<double, 1> u = Vector<double, 1>::Zero();
    filter.predict(u);

    Vector<double, 1> z;
    z << 1.0;
    filter.update(z);

    CHECK(filter.nees() >= 0.0);
    CHECK(std::isfinite(filter.nees()));
}
