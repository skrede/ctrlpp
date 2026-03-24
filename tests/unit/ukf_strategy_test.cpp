#include "ctrlpp/estimation/ukf.h"
#include "ctrlpp/estimation/ekf.h"
#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/sigma_points/julier_sigma_points.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

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

// Concept satisfaction static asserts
static_assert(ObserverPolicy<ukf<double, 2, 1, 1, ukf_linear_dynamics, ukf_position_measurement>>);
static_assert(CovarianceObserver<ukf<double, 2, 1, 1, ukf_linear_dynamics, ukf_position_measurement>>);

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("ukf with julier strategy")
{
    ukf_linear_dynamics dyn;
    ukf_position_measurement meas;

    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.01;
    Matrix<double, 1, 1> R;
    R << 1.0;
    Vector<double, 2> x0 = Vector<double, 2>::Zero();
    Matrix<double, 2, 2> P0 = Matrix<double, 2, 2>::Identity() * 10.0;

    ukf<double, 2, 1, 1, ukf_linear_dynamics, ukf_position_measurement, julier_sigma_points<double, 2>> filter(
        dyn, meas, ukf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0}, julier_options<double>{.kappa = 1.0});

    double true_pos = 0.0;
    double true_vel = 1.0;
    constexpr double dt = 0.1;

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
}

TEST_CASE("ukf with qr gain decomposition")
{
    ukf_linear_dynamics dyn;
    ukf_position_measurement meas;

    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.01;
    Matrix<double, 1, 1> R;
    R << 1.0;
    Vector<double, 2> x0 = Vector<double, 2>::Zero();
    Matrix<double, 2, 2> P0 = Matrix<double, 2, 2>::Identity() * 10.0;

    ukf filter(dyn, meas, ukf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0, .decomposition = gain_decomposition::qr});

    double true_pos = 0.0;
    double true_vel = 1.0;
    constexpr double dt = 0.1;

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
}

TEST_CASE("ukf shares dynamics_model with ekf")
{
    auto shared_dynamics = [](const Vector<double, 2>& x, const Vector<double, 1>& u) -> Vector<double, 2>
    {
        constexpr double dt = 0.1;
        Vector<double, 2> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    };

    auto shared_meas = [](const Vector<double, 2>& x) -> Vector<double, 1>
    {
        Vector<double, 1> z;
        z(0) = x(0);
        return z;
    };

    static_assert(dynamics_model<decltype(shared_dynamics), double, 2, 1>);

    // Same dynamics works with EKF
    ekf ekf_filter(shared_dynamics, shared_meas, ekf_config<double, 2, 1, 1>{});

    // Same dynamics works with UKF
    ukf ukf_filter(shared_dynamics, shared_meas, ukf_config<double, 2, 1, 1>{});

    Vector<double, 1> u = Vector<double, 1>::Zero();
    ekf_filter.predict(u);
    ukf_filter.predict(u);

    Vector<double, 1> z;
    z << 1.0;
    ekf_filter.update(z);
    ukf_filter.update(z);

    CHECK(std::isfinite(ekf_filter.state()(0)));
    CHECK(std::isfinite(ukf_filter.state()(0)));
}

TEST_CASE("ukf satisfies ObserverPolicy and CovarianceObserver")
{
    // Compile-time concept checks (also covered by static_assert at file scope)
    static_assert(ObserverPolicy<ukf<double, 2, 1, 1, ukf_linear_dynamics, ukf_position_measurement>>);
    static_assert(CovarianceObserver<ukf<double, 2, 1, 1, ukf_linear_dynamics, ukf_position_measurement>>);

    // Runtime check that the interface works
    ukf_linear_dynamics dyn;
    ukf_position_measurement meas;

    ukf filter(dyn, meas, ukf_config<double, 2, 1, 1>{});

    Vector<double, 1> u = Vector<double, 1>::Zero();
    filter.predict(u);

    Vector<double, 1> z;
    z << 1.0;
    filter.update(z);

    [[maybe_unused]] const auto& s = filter.state();
    [[maybe_unused]] const auto& P = filter.covariance();
    [[maybe_unused]] const auto& inn = filter.innovation();

    CHECK(std::isfinite(s(0)));
}
