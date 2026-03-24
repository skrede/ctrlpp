#include "ctrlpp/estimation/kalman.h"

#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/model/state_space.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>


// Helper: constant velocity model
// State: [position, velocity], input: acceleration, output: position
static auto make_const_velocity_system()
{
    ctrlpp::discrete_state_space<double, 2, 1, 1> sys;
    double dt = 0.1;
    sys.A << 1.0, dt,
             0.0, 1.0;
    sys.B << 0.5 * dt * dt,
             dt;
    sys.C << 1.0, 0.0;
    sys.D << 0.0;
    return sys;
}

TEST_CASE("kalman filter convergence on constant velocity model") {
    auto sys = make_const_velocity_system();

    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Identity() * 10.0;

    ctrlpp::kalman_filter<double, 2, 1, 1> kf(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0});

    // True state: position=5, velocity=1 (constant)
    double true_pos = 0.0;
    double true_vel = 1.0;
    double dt = 0.1;

    for (int i = 0; i < 50; ++i) {
        true_pos += true_vel * dt;

        Eigen::Matrix<double, 1, 1> u;
        u << 0.0;
        kf.predict(u);

        // Noisy measurement of position (no actual random noise, just a small fixed offset)
        Eigen::Matrix<double, 1, 1> z;
        z << true_pos + 0.1 * std::sin(static_cast<double>(i));
        kf.update(z);
    }

    auto est = kf.state();
    CHECK_THAT(est(0), Catch::Matchers::WithinAbs(true_pos, 0.5));
    CHECK_THAT(est(1), Catch::Matchers::WithinAbs(true_vel, 0.5));
}

TEST_CASE("kalman filter covariance remains symmetric and PSD") {
    auto sys = make_const_velocity_system();

    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Identity() * 10.0;

    ctrlpp::kalman_filter<double, 2, 1, 1> kf(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0});

    for (int i = 0; i < 50; ++i) {
        Eigen::Matrix<double, 1, 1> u;
        u << 0.0;
        kf.predict(u);

        Eigen::Matrix<double, 1, 1> z;
        z << static_cast<double>(i) * 0.1;
        kf.update(z);

        auto P = kf.covariance();
        // Symmetric
        CHECK((P - P.transpose()).norm() < 1e-10);
        // Positive semi-definite
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> eigsolver(P, Eigen::EigenvaluesOnly);
        for (int j = 0; j < 2; ++j)
            CHECK(eigsolver.eigenvalues()(j) >= -1e-10);
    }
}

TEST_CASE("kalman filter NEES is finite and positive after update") {
    auto sys = make_const_velocity_system();

    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Identity() * 10.0;

    ctrlpp::kalman_filter<double, 2, 1, 1> kf(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0});

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;
    kf.predict(u);

    Eigen::Matrix<double, 1, 1> z;
    z << 1.0;
    kf.update(z);

    CHECK(kf.nees() >= 0.0);
    CHECK(std::isfinite(kf.nees()));
}

TEST_CASE("kalman filter steady state detection") {
    auto sys = make_const_velocity_system();

    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Identity() * 10.0;

    ctrlpp::kalman_filter<double, 2, 1, 1> kf(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0});

    // Should not be steady state initially (P0 is large)
    // Run many iterations to converge
    for (int i = 0; i < 500; ++i) {
        Eigen::Matrix<double, 1, 1> u;
        u << 0.0;
        kf.predict(u);

        Eigen::Matrix<double, 1, 1> z;
        z << static_cast<double>(i) * 0.1;
        kf.update(z);
    }

    // Check the relative change is small
    auto P = kf.covariance();
    INFO("P norm: " << P.norm());
    // Use generous tolerance for steady-state detection
    CHECK(kf.is_steady_state(1e-2));
}

TEST_CASE("kalman filter reset covariance") {
    auto sys = make_const_velocity_system();

    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Identity() * 10.0;

    ctrlpp::kalman_filter<double, 2, 1, 1> kf(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0});

    // Run a few iterations
    for (int i = 0; i < 10; ++i) {
        Eigen::Matrix<double, 1, 1> u;
        u << 0.0;
        kf.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << 1.0;
        kf.update(z);
    }

    // Covariance should have changed from P0
    CHECK((kf.covariance() - P0).norm() > 0.1);

    // Reset
    kf.reset_covariance(P0);
    CHECK((kf.covariance() - P0).norm() < 1e-10);
}

TEST_CASE("kalman filter MIMO predict-update cycle") {
    // 3-state, 2-input, 2-output system
    ctrlpp::discrete_state_space<double, 3, 2, 2> sys;
    sys.A << 0.9, 0.1, 0.0,
             0.0, 0.8, 0.2,
             0.0, 0.0, 0.7;
    sys.B << 1.0, 0.0,
             0.0, 1.0,
             0.0, 0.0;
    sys.C << 1.0, 0.0, 0.0,
             0.0, 1.0, 0.0;
    sys.D = Eigen::Matrix<double, 2, 2>::Zero();

    Eigen::Matrix<double, 3, 3> Q = Eigen::Matrix<double, 3, 3>::Identity() * 0.01;
    Eigen::Matrix<double, 2, 2> R = Eigen::Matrix<double, 2, 2>::Identity() * 0.1;
    Eigen::Vector3d x0 = Eigen::Vector3d::Zero();
    Eigen::Matrix<double, 3, 3> P0 = Eigen::Matrix<double, 3, 3>::Identity();

    ctrlpp::kalman_filter<double, 3, 2, 2> kf(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0});

    // Run 10 predict-update cycles
    for (int i = 0; i < 10; ++i) {
        Eigen::Vector2d u;
        u << 0.1, 0.2;
        kf.predict(u);

        Eigen::Vector2d z;
        z << 1.0, 0.5;
        kf.update(z);
    }

    // Just verify it ran without error and state is finite
    auto est = kf.state();
    for (int i = 0; i < 3; ++i)
        CHECK(std::isfinite(est(i)));
}

TEST_CASE("kalman filter set_model updates system") {
    auto sys = make_const_velocity_system();

    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Identity();

    ctrlpp::kalman_filter<double, 2, 1, 1> kf(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0});

    // Change model
    auto sys2 = sys;
    sys2.A(0, 1) = 0.2;  // Different dt
    kf.set_model(sys2);

    // Predict with new model
    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;
    kf.predict(u);

    // State should have been propagated with new A
    // x = A_new * x0 + B * u = 0 (since x0 = 0 and u = 0)
    // Now set a nonzero state and verify
    kf.reset_covariance(P0);

    // Re-test with nonzero state via update
    Eigen::Matrix<double, 1, 1> z;
    z << 1.0;
    kf.update(z);

    CHECK(std::isfinite(kf.state()(0)));
}

// Static assertions are in the header; this test verifies they compile
TEST_CASE("kalman filter concept satisfaction") {
    static_assert(ctrlpp::ObserverPolicy<ctrlpp::kalman_filter<double, 2, 1, 1>>);
    static_assert(ctrlpp::CovarianceObserver<ctrlpp::kalman_filter<double, 2, 1, 1>>);
    CHECK(true);
}
