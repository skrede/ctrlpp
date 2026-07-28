#include "hardening_helpers.h"

#include "ctrlpp/mhe.h"
#include "ctrlpp/mpc/osqp_solver.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>

#include <cmath>
#include <cstddef>
#include <limits>
#include <random>

using Catch::Matchers::WithinAbs;

namespace
{

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr std::size_t NY = 1;
constexpr std::size_t N = 5;
constexpr double dt = 0.1;

struct linear_dynamics
{
    auto operator()(const ctrlpp::Vector<double, NX>& x, const ctrlpp::Vector<double, NU>& u) const -> ctrlpp::Vector<double, NX>
    {
        ctrlpp::Vector<double, NX> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    }

    auto jacobian_x(const ctrlpp::Vector<double, NX>&, const ctrlpp::Vector<double, NU>&) const -> ctrlpp::Matrix<double, NX, NX>
    {
        ctrlpp::Matrix<double, NX, NX> F;
        F << 1.0, dt, 0.0, 1.0;
        return F;
    }

    auto jacobian_u(const ctrlpp::Vector<double, NX>&, const ctrlpp::Vector<double, NU>&) const -> ctrlpp::Matrix<double, NX, NU>
    {
        ctrlpp::Matrix<double, NX, NU> G;
        G << 0.5 * dt * dt, dt;
        return G;
    }
};

struct position_measurement
{
    auto operator()(const ctrlpp::Vector<double, NX>& x) const -> ctrlpp::Vector<double, NY>
    {
        return (ctrlpp::Vector<double, NY>() << x(0)).finished();
    }

    auto jacobian(const ctrlpp::Vector<double, NX>&) const -> ctrlpp::Matrix<double, NY, NX>
    {
        return (ctrlpp::Matrix<double, NY, NX>() << 1.0, 0.0).finished();
    }
};

using MheType = ctrlpp::mhe<double, NX, NU, NY, N, ctrlpp::osqp_solver, linear_dynamics, position_measurement>;

}

// ── MHE hardening: negative ────────────────────────────────────────────────────

TEST_CASE("MHE with NaN in measurement noise", "[mhe][hardening][negative]")
{
    ctrlpp::mhe_config<double, NX, NU, NY, N> cfg;
    cfg.Q = ctrlpp::Matrix<double, NX, NX>::Identity() * 0.01;
    cfg.R = ctrlpp::Matrix<double, NY, NY>::Identity() * 0.1;

    auto estimator = ctrlpp::test::constructed(MheType::create(linear_dynamics{}, position_measurement{}, cfg));

    ctrlpp::Vector<double, NU> u = ctrlpp::Vector<double, NU>::Zero();

    // Warm up with valid measurements
    for (std::size_t i = 0; i < N + 1; ++i) {
        estimator.predict(u);
        ctrlpp::Vector<double, NY> z;
        z << 0.1 * static_cast<double>(i);
        estimator.update(z);
    }

    // Inject NaN measurement
    estimator.predict(u);
    ctrlpp::Vector<double, NY> z_nan;
    z_nan << std::numeric_limits<double>::quiet_NaN();
    estimator.update(z_nan);

    auto x_hat = estimator.state();
    // NaN should propagate -- state may become NaN but should not crash
    (void)x_hat;
}

TEST_CASE("MHE with inconsistent measurements", "[mhe][hardening][negative]")
{
    ctrlpp::mhe_config<double, NX, NU, NY, N> cfg;
    cfg.Q = ctrlpp::Matrix<double, NX, NX>::Identity() * 0.01;
    cfg.R = ctrlpp::Matrix<double, NY, NY>::Identity() * 0.1;

    auto estimator = ctrlpp::test::constructed(MheType::create(linear_dynamics{}, position_measurement{}, cfg));

    ctrlpp::Vector<double, NU> u = ctrlpp::Vector<double, NU>::Zero();

    // Feed wildly inconsistent measurements (jumping between +100 and -100)
    for (std::size_t i = 0; i < 2 * N; ++i) {
        estimator.predict(u);
        ctrlpp::Vector<double, NY> z;
        z << ((i % 2 == 0) ? 100.0 : -100.0);
        estimator.update(z);
    }

    auto x_hat = estimator.state();
    REQUIRE(std::isfinite(x_hat(0)));
    REQUIRE(std::isfinite(x_hat(1)));
}

// ── MHE hardening: precision ───────────────────────────────────────────────────

TEST_CASE("MHE linear system matches Kalman-like estimate", "[mhe][hardening][precision]")
{
    ctrlpp::mhe_config<double, NX, NU, NY, N> cfg;
    cfg.Q = ctrlpp::Matrix<double, NX, NX>::Identity() * 0.01;
    cfg.R = ctrlpp::Matrix<double, NY, NY>::Identity() * 0.1;
    cfg.P0 = ctrlpp::Matrix<double, NX, NX>::Identity() * 10.0;

    auto estimator = ctrlpp::test::constructed(MheType::create(linear_dynamics{}, position_measurement{}, cfg));

    // True state: position 0, velocity 1 (constant velocity)
    ctrlpp::Vector<double, NX> x_true;
    x_true << 0.0, 1.0;
    ctrlpp::Vector<double, NU> u = ctrlpp::Vector<double, NU>::Zero();

    std::mt19937 gen(42);
    std::normal_distribution<double> noise(0.0, 0.1);

    for (std::size_t i = 0; i < 3 * N; ++i) {
        x_true = linear_dynamics{}(x_true, u);
        ctrlpp::Vector<double, NY> z;
        z << x_true(0) + noise(gen);

        estimator.predict(u);
        estimator.update(z);
    }

    auto x_hat = estimator.state();
    // Position estimate should be close to truth
    REQUIRE_THAT(x_hat(0), WithinAbs(x_true(0), 0.5));
}

// ── MHE hardening: convergence ─────────────────────────────────────────────────

TEST_CASE("MHE state estimate converges to truth", "[mhe][hardening][convergence]")
{
    ctrlpp::mhe_config<double, NX, NU, NY, N> cfg;
    cfg.Q = ctrlpp::Matrix<double, NX, NX>::Identity() * 0.01;
    cfg.R = ctrlpp::Matrix<double, NY, NY>::Identity() * 0.1;
    cfg.P0 = ctrlpp::Matrix<double, NX, NX>::Identity() * 100.0;
    cfg.x0 = ctrlpp::Vector<double, NX>::Zero(); // start far from truth

    auto estimator = ctrlpp::test::constructed(MheType::create(linear_dynamics{}, position_measurement{}, cfg));

    ctrlpp::Vector<double, NX> x_true;
    x_true << 5.0, 0.5;
    ctrlpp::Vector<double, NU> u = ctrlpp::Vector<double, NU>::Zero();

    std::mt19937 gen(99);
    std::normal_distribution<double> noise(0.0, 0.1);

    double prev_error = 1e10;
    int improvements = 0;

    for (std::size_t i = 0; i < 5 * N; ++i) {
        x_true = linear_dynamics{}(x_true, u);
        ctrlpp::Vector<double, NY> z;
        z << x_true(0) + noise(gen);

        estimator.predict(u);
        estimator.update(z);

        double error = (estimator.state() - x_true).norm();
        if (error < prev_error) {
            ++improvements;
        }
        prev_error = error;
    }

    // Over time, error should generally decrease
    CHECK(improvements > static_cast<int>(2 * N));
    CHECK(prev_error < 5.0);
}

// ── MHE hardening: robustness ──────────────────────────────────────────────────

TEST_CASE("MHE with ill-conditioned process noise", "[mhe][hardening][robustness]")
{
    ctrlpp::mhe_config<double, NX, NU, NY, N> cfg;
    auto Q_ill = ctrlpp::test::ill_conditioned_2x2<double>(1e10);
    cfg.Q = Q_ill;
    cfg.R = ctrlpp::Matrix<double, NY, NY>::Identity() * 0.1;

    auto estimator = ctrlpp::test::constructed(MheType::create(linear_dynamics{}, position_measurement{}, cfg));

    ctrlpp::Vector<double, NU> u = ctrlpp::Vector<double, NU>::Zero();

    for (std::size_t i = 0; i < 2 * N; ++i) {
        estimator.predict(u);
        ctrlpp::Vector<double, NY> z;
        z << static_cast<double>(i) * 0.1;
        estimator.update(z);
    }

    auto x_hat = estimator.state();
    REQUIRE(std::isfinite(x_hat(0)));
    REQUIRE(std::isfinite(x_hat(1)));
}
