#include "hardening_helpers.h"
#include "ctrlpp/estimation/ekf.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <numbers>

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

struct position_measurement
{
    auto operator()(const Vector<double, 2>& x) const -> Vector<double, 1>
    {
        return (Vector<double, 1>() << x(0)).finished();
    }

    auto jacobian(const Vector<double, 2>& /*x*/) const -> Matrix<double, 1, 2>
    {
        return (Matrix<double, 1, 2>() << 1.0, 0.0).finished();
    }
};

// ---------------------------------------------------------------------------
// Test dynamics/measurement: pendulum with analytical Jacobians
// ---------------------------------------------------------------------------
struct pendulum_dynamics
{
    static constexpr double g = 9.81;
    static constexpr double l = 1.0;
    static constexpr double b = 0.1;
    static constexpr double m = 1.0;
    static constexpr double dt = 0.01;

    auto operator()(const Vector<double, 2>& x, const Vector<double, 1>& u) const -> Vector<double, 2>
    {
        double theta = x(0), omega = x(1), tau = u(0);
        Vector<double, 2> x_next;
        x_next(0) = theta + omega * dt;
        x_next(1) = omega + (-g / l * std::sin(theta) - b * omega + tau / (m * l * l)) * dt;
        return x_next;
    }

    auto jacobian_x(const Vector<double, 2>& x, const Vector<double, 1>& /*u*/) const -> Matrix<double, 2, 2>
    {
        Matrix<double, 2, 2> F;
        F << 1.0, dt, -g / l * std::cos(x(0)) * dt, 1.0 - b * dt;
        return F;
    }

    auto jacobian_u(const Vector<double, 2>& /*x*/, const Vector<double, 1>& /*u*/) const -> Matrix<double, 2, 1>
    {
        return (Matrix<double, 2, 1>() << 0.0, dt / (m * l * l)).finished();
    }
};

struct angle_measurement
{
    auto operator()(const Vector<double, 2>& x) const -> Vector<double, 1>
    {
        return (Vector<double, 1>() << x(0)).finished();
    }

    auto jacobian(const Vector<double, 2>& /*x*/) const -> Matrix<double, 1, 2>
    {
        return (Matrix<double, 1, 2>() << 1.0, 0.0).finished();
    }
};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("ekf with analytical Jacobians converges on linear system")
{
    linear_dynamics dyn;
    position_measurement meas;

    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.01;
    Matrix<double, 1, 1> R;
    R << 1.0;
    Vector<double, 2> x0 = Vector<double, 2>::Zero();
    Matrix<double, 2, 2> P0 = Matrix<double, 2, 2>::Identity() * 10.0;

    auto filter = ctrlpp::test::constructed(ekf<double, 2, 1, 1, decltype(dyn), decltype(meas)>::create(dyn, meas, ekf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0}));

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
        REQUIRE(filter.update(z).has_value());
    }

    auto est = filter.state();
    CHECK_THAT(est(0), Catch::Matchers::WithinAbs(true_pos, 0.5));
    CHECK_THAT(est(1), Catch::Matchers::WithinAbs(true_vel, 0.5));
}

TEST_CASE("ekf with numerical Jacobians converges on linear system")
{
    auto dyn = [](const Vector<double, 2>& x, const Vector<double, 1>& u) -> Vector<double, 2>
    {
        constexpr double dt = 0.1;
        Vector<double, 2> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    };

    auto meas = [](const Vector<double, 2>& x) -> Vector<double, 1>
    {
        Vector<double, 1> z;
        z(0) = x(0);
        return z;
    };

    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.01;
    Matrix<double, 1, 1> R;
    R << 1.0;
    Vector<double, 2> x0 = Vector<double, 2>::Zero();
    Matrix<double, 2, 2> P0 = Matrix<double, 2, 2>::Identity() * 10.0;

    auto filter = ctrlpp::test::constructed(ekf<double, 2, 1, 1, decltype(dyn), decltype(meas)>::create(dyn, meas, ekf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0}));

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
        REQUIRE(filter.update(z).has_value());
    }

    auto est = filter.state();
    CHECK_THAT(est(0), Catch::Matchers::WithinAbs(true_pos, 0.5));
    CHECK_THAT(est(1), Catch::Matchers::WithinAbs(true_vel, 0.5));

    // Both the dynamics and measurement above are linear, so the exact
    // analytic Jacobians F, G, H used by this same recursion are known
    // closed forms; the standard linear Kalman recursion built from those
    // closed forms is therefore an external, independent oracle for this
    // filter's central-difference numerical Jacobian path. Two independent
    // implementations of that recursion (numpy/scipy and Octave, see
    // validation/oracles/linear_kalman/oracle_kf.{py,m}) agree on
    // the state after 50 predict/update steps to within a few ULP
    // (observed max |disagreement| ~ 1e-15, i.e. within a handful of
    // std::numeric_limits<double>::epsilon()), which is the expected
    // agreement for two runs of the same associativity-sensitive
    // floating-point recursion; the shared value below is that
    // dual-oracle-agreed golden, not ctrlpp's own output.
    //
    // The central-difference stencil (detail/numerical_diff.h) applies a
    // cbrt(eps) step; for an exactly linear map the truncation term
    // vanishes and only the round-off floor eps^(2/3) remains per Jacobian
    // entry (same accuracy floor derived and asserted analytically in
    // numerical_diff_test.cpp). That per-entry floor accumulates additively
    // and non-adversarially across the state dimension and the predict/
    // update recursion length, bounded by the problem's own state-magnitude
    // scale (the true position grows to n_steps * dt * true_vel = 5 over
    // the run):
    constexpr int n_steps = 50;
    constexpr int nx = 2;
    const double eps = std::numeric_limits<double>::epsilon();
    const double jacobian_step_tol = static_cast<double>(n_steps) * static_cast<double>(nx) * std::pow(eps, 2.0 / 3.0) * (n_steps * dt * true_vel);

    const double oracle_pos = 4.9849480216340529;
    const double oracle_vel = 0.99088209456807685;
    CHECK_THAT(est(0), Catch::Matchers::WithinAbs(oracle_pos, jacobian_step_tol));
    CHECK_THAT(est(1), Catch::Matchers::WithinAbs(oracle_vel, jacobian_step_tol));
}

TEST_CASE("ekf nonlinear pendulum tracking")
{
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

    auto filter = ctrlpp::test::constructed(ekf<double, 2, 1, 1, decltype(dyn), decltype(meas)>::create(dyn, meas, ekf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0}));

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
        REQUIRE(filter.update(z).has_value());
    }

    auto est = filter.state();
    CHECK_THAT(est(0), Catch::Matchers::WithinAbs(x_true(0), 0.3));
    CHECK_THAT(est(1), Catch::Matchers::WithinAbs(x_true(1), 0.5));
}

TEST_CASE("ekf with shared dynamics_model lambda compiles and runs")
{
    auto shared_dynamics = [](const Vector<double, 2>& x, const Vector<double, 1>& u) -> Vector<double, 2>
    {
        constexpr double dt = 0.1;
        Vector<double, 2> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    };

    auto meas = [](const Vector<double, 2>& x) -> Vector<double, 1>
    {
        Vector<double, 1> z;
        z(0) = x(0);
        return z;
    };

    static_assert(dynamics_model<decltype(shared_dynamics), double, 2, 1>);

    auto filter = ctrlpp::test::constructed(ekf<double, 2, 1, 1, decltype(shared_dynamics), decltype(meas)>::create(shared_dynamics, meas, ekf_config<double, 2, 1, 1>{}));

    Vector<double, 1> u = Vector<double, 1>::Zero();
    filter.predict(u);

    Vector<double, 1> z;
    z << 1.0;
    REQUIRE(filter.update(z).has_value());

    CHECK(std::isfinite(filter.state()(0)));
}
