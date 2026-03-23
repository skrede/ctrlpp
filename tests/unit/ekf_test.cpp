#include "ctrlpp/ekf.h"

#include "ctrlpp/observer_policy.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>
#include <numbers>

using namespace ctrlpp;

// ---------------------------------------------------------------------------
// Test dynamics: linear constant-velocity model with analytical Jacobians
// ---------------------------------------------------------------------------
struct linear_dynamics {
    double dt = 0.1;

    auto operator()(const Vector<double, 2>& x,
                    const Vector<double, 1>& u) const -> Vector<double, 2>
    {
        Vector<double, 2> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    }

    auto jacobian_x(const Vector<double, 2>& /*x*/,
                    const Vector<double, 1>& /*u*/) const -> Matrix<double, 2, 2>
    {
        Matrix<double, 2, 2> F;
        F << 1.0, dt,
             0.0, 1.0;
        return F;
    }

    auto jacobian_u(const Vector<double, 2>& /*x*/,
                    const Vector<double, 1>& /*u*/) const -> Matrix<double, 2, 1>
    {
        Matrix<double, 2, 1> G;
        G << 0.5 * dt * dt,
             dt;
        return G;
    }
};

// ---------------------------------------------------------------------------
// Test measurement: position observation with analytical Jacobian
// ---------------------------------------------------------------------------
struct position_measurement {
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
// Test dynamics/measurement: pendulum with analytical Jacobians
// ---------------------------------------------------------------------------
struct pendulum_dynamics {
    static constexpr double g = 9.81;
    static constexpr double l = 1.0;
    static constexpr double b = 0.1;
    static constexpr double m = 1.0;
    static constexpr double dt = 0.01;

    auto operator()(const Vector<double, 2>& x,
                    const Vector<double, 1>& u) const -> Vector<double, 2>
    {
        Vector<double, 2> x_next;
        double theta = x(0);
        double omega = x(1);
        double tau = u(0);
        x_next(0) = theta + omega * dt;
        x_next(1) = omega + (-g / l * std::sin(theta) - b * omega + tau / (m * l * l)) * dt;
        return x_next;
    }

    auto jacobian_x(const Vector<double, 2>& x,
                    const Vector<double, 1>& /*u*/) const -> Matrix<double, 2, 2>
    {
        Matrix<double, 2, 2> F;
        F << 1.0, dt,
             -g / l * std::cos(x(0)) * dt, 1.0 - b * dt;
        return F;
    }

    auto jacobian_u(const Vector<double, 2>& /*x*/,
                    const Vector<double, 1>& /*u*/) const -> Matrix<double, 2, 1>
    {
        Matrix<double, 2, 1> G;
        G << 0.0,
             dt / (m * l * l);
        return G;
    }
};

struct angle_measurement {
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
// Concept satisfaction static asserts
// ---------------------------------------------------------------------------
static_assert(ObserverPolicy<ekf<double, 2, 1, 1, linear_dynamics, position_measurement>>);
static_assert(CovarianceObserver<ekf<double, 2, 1, 1, linear_dynamics, position_measurement>>);

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

TEST_CASE("ekf with analytical Jacobians converges on linear system") {
    linear_dynamics dyn;
    position_measurement meas;

    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.01;
    Matrix<double, 1, 1> R;
    R << 1.0;
    Vector<double, 2> x0 = Vector<double, 2>::Zero();
    Matrix<double, 2, 2> P0 = Matrix<double, 2, 2>::Identity() * 10.0;

    ekf filter(dyn, meas, ekf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0});

    double true_pos = 0.0;
    double true_vel = 1.0;
    constexpr double dt = 0.1;

    for (int i = 0; i < 50; ++i) {
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

TEST_CASE("ekf with numerical Jacobians converges on linear system") {
    auto dyn = [](const Vector<double, 2>& x,
                  const Vector<double, 1>& u) -> Vector<double, 2> {
        constexpr double dt = 0.1;
        Vector<double, 2> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    };

    auto meas = [](const Vector<double, 2>& x) -> Vector<double, 1> {
        Vector<double, 1> z;
        z(0) = x(0);
        return z;
    };

    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.01;
    Matrix<double, 1, 1> R;
    R << 1.0;
    Vector<double, 2> x0 = Vector<double, 2>::Zero();
    Matrix<double, 2, 2> P0 = Matrix<double, 2, 2>::Identity() * 10.0;

    ekf filter(dyn, meas, ekf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0});

    double true_pos = 0.0;
    double true_vel = 1.0;
    constexpr double dt = 0.1;

    for (int i = 0; i < 50; ++i) {
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

TEST_CASE("ekf nonlinear pendulum tracking") {
    pendulum_dynamics dyn;
    angle_measurement meas;

    constexpr double dt = pendulum_dynamics::dt;
    constexpr double g = pendulum_dynamics::g;
    constexpr double l = pendulum_dynamics::l;
    constexpr double b = pendulum_dynamics::b;

    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.001;
    Matrix<double, 1, 1> R;
    R << 0.01;
    Vector<double, 2> x0;
    x0 << 0.1, 0.0;
    Matrix<double, 2, 2> P0 = Matrix<double, 2, 2>::Identity() * 1.0;

    ekf filter(dyn, meas, ekf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0});

    Vector<double, 2> x_true;
    x_true << std::numbers::pi / 4.0, 0.0;

    for (int i = 0; i < 200; ++i) {
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

TEST_CASE("ekf covariance stays symmetric and PSD over 100+ cycles") {
    linear_dynamics dyn;
    position_measurement meas;

    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.01;
    Matrix<double, 1, 1> R;
    R << 1.0;
    Vector<double, 2> x0 = Vector<double, 2>::Zero();
    Matrix<double, 2, 2> P0 = Matrix<double, 2, 2>::Identity() * 10.0;

    ekf filter(dyn, meas, ekf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0});

    for (int i = 0; i < 150; ++i) {
        Vector<double, 1> u = Vector<double, 1>::Zero();
        filter.predict(u);

        Vector<double, 1> z;
        z << static_cast<double>(i) * 0.1;
        filter.update(z);

        auto P = filter.covariance();
        CHECK((P - P.transpose()).norm() < 1e-10);
        Eigen::SelfAdjointEigenSolver<Matrix<double, 2, 2>> eigsolver(P, Eigen::EigenvaluesOnly);
        for (int j = 0; j < 2; ++j)
            CHECK(eigsolver.eigenvalues()(j) >= -1e-10);
    }
}

TEST_CASE("ekf NEES is finite and positive") {
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

TEST_CASE("ekf with shared dynamics_model lambda compiles and runs") {
    auto shared_dynamics = [](const Vector<double, 2>& x,
                              const Vector<double, 1>& u) -> Vector<double, 2> {
        constexpr double dt = 0.1;
        Vector<double, 2> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    };

    auto meas = [](const Vector<double, 2>& x) -> Vector<double, 1> {
        Vector<double, 1> z;
        z(0) = x(0);
        return z;
    };

    static_assert(dynamics_model<decltype(shared_dynamics), double, 2, 1>);

    ekf filter(shared_dynamics, meas, ekf_config<double, 2, 1, 1>{});

    Vector<double, 1> u = Vector<double, 1>::Zero();
    filter.predict(u);

    Vector<double, 1> z;
    z << 1.0;
    filter.update(z);

    CHECK(std::isfinite(filter.state()(0)));
}
