#include "ctrlpp/estimation/ukf.h"
#include "ctrlpp/estimation/ekf.h"
#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/sigma_points/merwe_sigma_points.h"
#include "ctrlpp/estimation/sigma_points/julier_sigma_points.h"
#include "ctrlpp/estimation/sigma_points/sigma_point_strategy.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>
#include <numbers>
#include <numeric>

using namespace ctrlpp;

// ---------------------------------------------------------------------------
// Test dynamics: linear constant-velocity (double integrator)
// ---------------------------------------------------------------------------
struct ukf_linear_dynamics {
    double dt = 0.1;

    auto operator()(const Vector<double, 2>& x,
                    const Vector<double, 1>& u) const -> Vector<double, 2>
    {
        Vector<double, 2> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    }
};

// ---------------------------------------------------------------------------
// Test measurement: position observation
// ---------------------------------------------------------------------------
struct ukf_position_measurement {
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
struct ukf_pendulum_dynamics {
    static constexpr double g = 9.81;
    static constexpr double l = 1.0;
    static constexpr double b = 0.1;
    static constexpr double m = 1.0;
    static constexpr double dt = 0.01;

    auto operator()(const Vector<double, 2>& x,
                    const Vector<double, 1>& u) const -> Vector<double, 2>
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

struct ukf_angle_measurement {
    auto operator()(const Vector<double, 2>& x) const -> Vector<double, 1>
    {
        Vector<double, 1> z;
        z(0) = x(0);
        return z;
    }
};

// ---------------------------------------------------------------------------
// Concept satisfaction static asserts
// ---------------------------------------------------------------------------
static_assert(sigma_point_strategy<merwe_sigma_points<double, 2>, double, 2>);
static_assert(sigma_point_strategy<julier_sigma_points<double, 2>, double, 2>);

static_assert(ObserverPolicy<ukf<double, 2, 1, 1, ukf_linear_dynamics, ukf_position_measurement>>);
static_assert(CovarianceObserver<ukf<double, 2, 1, 1, ukf_linear_dynamics, ukf_position_measurement>>);

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("merwe sigma points generate correct weights") {
    constexpr std::size_t NX = 2;
    merwe_sigma_points<double, NX> sp{merwe_options<double>{.alpha = 1e-3, .beta = 2.0, .kappa = 0.0}};

    Vector<double, NX> x;
    x << 1.0, 2.0;
    Matrix<double, NX, NX> P = Matrix<double, NX, NX>::Identity();

    auto result = sp.generate(x, P);

    // Wm should sum to 1.0
    double wm_sum = 0.0;
    for (auto w : result.Wm) wm_sum += w;
    CHECK_THAT(wm_sum, Catch::Matchers::WithinAbs(1.0, 1e-10));

    // Wc sum: Wc_0 = Wm_0 + (1 - alpha^2 + beta), so Wc sums to Wm_sum + (1 - alpha^2 + beta)
    double wc_sum = 0.0;
    for (auto w : result.Wc) wc_sum += w;
    double expected_wc_sum = 1.0 + (1.0 - 1e-3 * 1e-3 + 2.0);
    CHECK_THAT(wc_sum, Catch::Matchers::WithinAbs(expected_wc_sum, 1e-10));

    // Verify Wm_0 and Wc_0 formulas
    double alpha = 1e-3;
    double beta = 2.0;
    double kappa = 0.0;
    double n = static_cast<double>(NX);
    double lambda = alpha * alpha * (n + kappa) - n;
    double expected_wm0 = lambda / (n + lambda);
    double expected_wc0 = expected_wm0 + (1.0 - alpha * alpha + beta);

    CHECK_THAT(result.Wm[0], Catch::Matchers::WithinAbs(expected_wm0, 1e-12));
    CHECK_THAT(result.Wc[0], Catch::Matchers::WithinAbs(expected_wc0, 1e-12));

    // Verify Wm_i = Wc_i = 1/(2*(n+lambda))
    double expected_wi = 1.0 / (2.0 * (n + lambda));
    for (std::size_t i = 1; i < merwe_sigma_points<double, NX>::num_points; ++i) {
        CHECK_THAT(result.Wm[i], Catch::Matchers::WithinAbs(expected_wi, 1e-12));
        CHECK_THAT(result.Wc[i], Catch::Matchers::WithinAbs(expected_wi, 1e-12));
    }
}

TEST_CASE("merwe sigma points capture mean and covariance") {
    constexpr std::size_t NX = 2;
    merwe_sigma_points<double, NX> sp{merwe_options<double>{.alpha = 1e-1, .beta = 2.0, .kappa = 0.0}};

    Vector<double, NX> x;
    x << 3.0, -1.0;
    Matrix<double, NX, NX> P;
    P << 4.0, 1.0,
         1.0, 2.0;

    auto result = sp.generate(x, P);

    // Weighted mean should reconstruct x
    Vector<double, NX> x_recon = Vector<double, NX>::Zero();
    for (std::size_t i = 0; i < merwe_sigma_points<double, NX>::num_points; ++i) {
        x_recon += result.Wm[i] * result.points[i];
    }
    CHECK_THAT(x_recon(0), Catch::Matchers::WithinAbs(x(0), 1e-10));
    CHECK_THAT(x_recon(1), Catch::Matchers::WithinAbs(x(1), 1e-10));

    // Weighted covariance should reconstruct P
    Matrix<double, NX, NX> P_recon = Matrix<double, NX, NX>::Zero();
    for (std::size_t i = 0; i < merwe_sigma_points<double, NX>::num_points; ++i) {
        auto diff = (result.points[i] - x_recon).eval();
        P_recon += result.Wc[i] * diff * diff.transpose();
    }
    CHECK_THAT(P_recon(0, 0), Catch::Matchers::WithinAbs(P(0, 0), 1e-8));
    CHECK_THAT(P_recon(0, 1), Catch::Matchers::WithinAbs(P(0, 1), 1e-8));
    CHECK_THAT(P_recon(1, 0), Catch::Matchers::WithinAbs(P(1, 0), 1e-8));
    CHECK_THAT(P_recon(1, 1), Catch::Matchers::WithinAbs(P(1, 1), 1e-8));
}

TEST_CASE("julier sigma points satisfy concept and generate") {
    constexpr std::size_t NX = 2;
    static_assert(sigma_point_strategy<julier_sigma_points<double, NX>, double, NX>);

    julier_sigma_points<double, NX> sp{julier_options<double>{.kappa = 0.0}};

    Vector<double, NX> x;
    x << 1.0, 2.0;
    Matrix<double, NX, NX> P = Matrix<double, NX, NX>::Identity();

    auto result = sp.generate(x, P);

    // Wm should sum to 1.0
    double wm_sum = 0.0;
    for (auto w : result.Wm) wm_sum += w;
    CHECK_THAT(wm_sum, Catch::Matchers::WithinAbs(1.0, 1e-10));

    // Wc should sum to 1.0
    double wc_sum = 0.0;
    for (auto w : result.Wc) wc_sum += w;
    CHECK_THAT(wc_sum, Catch::Matchers::WithinAbs(1.0, 1e-10));

    // Verify Wm_0 = Wc_0 = kappa/(NX + kappa)
    double n = static_cast<double>(NX);
    double kappa = 0.0;
    double expected_w0 = kappa / (n + kappa);
    CHECK_THAT(result.Wm[0], Catch::Matchers::WithinAbs(expected_w0, 1e-12));
    CHECK_THAT(result.Wc[0], Catch::Matchers::WithinAbs(expected_w0, 1e-12));

    // Verify Wm_i = Wc_i = 1/(2*(n+kappa))
    double expected_wi = 1.0 / (2.0 * (n + kappa));
    for (std::size_t i = 1; i < julier_sigma_points<double, NX>::num_points; ++i) {
        CHECK_THAT(result.Wm[i], Catch::Matchers::WithinAbs(expected_wi, 1e-12));
        CHECK_THAT(result.Wc[i], Catch::Matchers::WithinAbs(expected_wi, 1e-12));
    }
}

TEST_CASE("ukf tracks linear system") {
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

    // Covariance should have converged (trace decreased)
    CHECK(filter.covariance().trace() < initial_trace);
}

TEST_CASE("ukf tracks nonlinear system") {
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

TEST_CASE("ukf covariance remains PSD") {
    ukf_linear_dynamics dyn;
    ukf_position_measurement meas;

    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.01;
    Matrix<double, 1, 1> R;
    R << 1.0;
    Vector<double, 2> x0 = Vector<double, 2>::Zero();
    Matrix<double, 2, 2> P0 = Matrix<double, 2, 2>::Identity() * 10.0;

    ukf filter(dyn, meas, ukf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0});

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

TEST_CASE("ukf with julier strategy") {
    ukf_linear_dynamics dyn;
    ukf_position_measurement meas;

    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.01;
    Matrix<double, 1, 1> R;
    R << 1.0;
    Vector<double, 2> x0 = Vector<double, 2>::Zero();
    Matrix<double, 2, 2> P0 = Matrix<double, 2, 2>::Identity() * 10.0;

    ukf<double, 2, 1, 1, ukf_linear_dynamics, ukf_position_measurement,
        julier_sigma_points<double, 2>> filter(
        dyn, meas,
        ukf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0},
        julier_options<double>{.kappa = 1.0});

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

TEST_CASE("ukf with qr gain decomposition") {
    ukf_linear_dynamics dyn;
    ukf_position_measurement meas;

    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.01;
    Matrix<double, 1, 1> R;
    R << 1.0;
    Vector<double, 2> x0 = Vector<double, 2>::Zero();
    Matrix<double, 2, 2> P0 = Matrix<double, 2, 2>::Identity() * 10.0;

    ukf filter(dyn, meas, ukf_config<double, 2, 1, 1>{
        .Q = Q, .R = R, .x0 = x0, .P0 = P0,
        .decomposition = gain_decomposition::qr});

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

TEST_CASE("ukf shares dynamics_model with ekf") {
    auto shared_dynamics = [](const Vector<double, 2>& x,
                              const Vector<double, 1>& u) -> Vector<double, 2> {
        constexpr double dt = 0.1;
        Vector<double, 2> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    };

    auto shared_meas = [](const Vector<double, 2>& x) -> Vector<double, 1> {
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

TEST_CASE("ukf satisfies ObserverPolicy and CovarianceObserver") {
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
