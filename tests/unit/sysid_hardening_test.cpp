#include "hardening_helpers.h"

#include "ctrlpp/sysid/rls.h"
#include "ctrlpp/sysid/batch_arx.h"
#include "ctrlpp/sysid/moesp.h"
#include "ctrlpp/sysid/recursive_arx.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <random>

using Catch::Matchers::WithinAbs;

// ── RLS hardening ──────────────────────────────────────────────────────────────

TEST_CASE("RLS with singular covariance P0", "[rls][hardening][negative]")
{
    ctrlpp::rls_config<double, 2> cfg;
    cfg.P0 = ctrlpp::Matrix<double, 2, 2>::Zero();
    ctrlpp::rls<double, 2> estimator(cfg);

    Eigen::Vector2d phi;
    phi << 1.0, 0.5;
    estimator.update(1.0, phi);

    auto theta = estimator.parameters();
    REQUIRE(std::isfinite(theta(0)));
    REQUIRE(std::isfinite(theta(1)));
}

TEST_CASE("RLS with zero forgetting factor", "[rls][hardening][negative]")
{
    ctrlpp::rls_config<double, 2> cfg;
    cfg.lambda = 0.0;
    ctrlpp::rls<double, 2> estimator(cfg);

    Eigen::Vector2d phi;
    phi << 1.0, 0.5;

    // Zero lambda causes division by lambda in covariance update
    // Should not crash; result may be NaN/Inf but must not segfault
    estimator.update(1.0, phi);
    auto theta = estimator.parameters();
    // No crash is the success criterion
    (void)theta;
}

TEST_CASE("RLS with NaN regressor", "[rls][hardening][negative]")
{
    ctrlpp::rls<double, 2> estimator;

    auto phi = ctrlpp::test::nan_vector<double, 2>();
    estimator.update(1.0, phi);

    auto theta = estimator.parameters();
    // NaN regressor causes degenerate denominator -- update is skipped,
    // parameters remain at their initial value (zero)
    CHECK(theta(0) == 0.0);
    CHECK(theta(1) == 0.0);
}

TEST_CASE("RLS identifies known first-order system", "[rls][hardening][convergence]")
{
    constexpr std::size_t NP = 2;
    ctrlpp::rls<double, NP> estimator;

    std::mt19937 gen(42);
    std::normal_distribution<double> noise(0.0, 0.001);
    std::uniform_real_distribution<double> input(-1.0, 1.0);

    Eigen::Vector2d true_theta;
    true_theta << 0.8, 0.5; // y = 0.8*x1 + 0.5*x2

    for (int i = 0; i < 500; ++i) {
        Eigen::Vector2d phi;
        phi << input(gen), input(gen);
        double y = true_theta.dot(phi) + noise(gen);
        estimator.update(y, phi);
    }

    auto theta_hat = estimator.parameters();
    REQUIRE_THAT(theta_hat(0), WithinAbs(0.8, 0.05));
    REQUIRE_THAT(theta_hat(1), WithinAbs(0.5, 0.05));
}

TEST_CASE("RLS with ill-conditioned regressor", "[rls][hardening][robustness]")
{
    ctrlpp::rls<double, 2> estimator;

    Eigen::Vector2d true_theta;
    true_theta << 1.0, 2.0;

    std::mt19937 gen(77);
    std::normal_distribution<double> noise(0.0, 0.01);

    for (int i = 0; i < 200; ++i) {
        // Ill-conditioned: second regressor element is 1e-10 * first
        double x = static_cast<double>(i) * 0.01;
        Eigen::Vector2d phi;
        phi << x, x * 1e-10;
        double y = true_theta.dot(phi) + noise(gen);
        estimator.update(y, phi);
    }

    auto theta = estimator.parameters();
    REQUIRE(std::isfinite(theta(0)));
    REQUIRE(std::isfinite(theta(1)));
}

// ── ARX hardening ──────────────────────────────────────────────────────────────

TEST_CASE("Batch ARX with insufficient data", "[arx][hardening][negative]")
{
    // Only 2 data points for NA=2, NB=1 model -- need at least max(NA,NB)+1
    Eigen::RowVectorXd Y(3);
    Y << 1.0, 2.0, 3.0;
    Eigen::RowVectorXd U(3);
    U << 0.5, 1.0, 1.5;

    // Should complete without crashing; results may be poor
    auto result = ctrlpp::batch_arx<2, 1>(Y, U);
    REQUIRE(std::isfinite(result.system.A(0, 0)));
}

TEST_CASE("Batch ARX identifies known AR/X coefficients", "[arx][hardening][convergence]")
{
    // Generate data from known ARX(1,1): y(t) = 0.7*y(t-1) + 0.3*u(t-1)
    constexpr int N = 500;
    Eigen::RowVectorXd Y(N);
    Eigen::RowVectorXd U(N);

    std::mt19937 gen(42);
    std::normal_distribution<double> noise(0.0, 0.001);
    std::uniform_real_distribution<double> input(-1.0, 1.0);

    Y(0) = 0.0;
    U(0) = input(gen);
    for (int t = 1; t < N; ++t) {
        U(t) = input(gen);
        Y(t) = 0.7 * Y(t - 1) + 0.3 * U(t - 1) + noise(gen);
    }

    auto result = ctrlpp::batch_arx<1, 1>(Y, U);
    // NRMSE close to 0 indicates good fit (norm_error / norm_centered)
    REQUIRE(result.metrics.nrmse < 0.1);
}

TEST_CASE("Batch ARX with rank-deficient regressors", "[arx][hardening][negative]")
{
    // Constant input and output -- rank-deficient regressor matrix
    constexpr int N = 50;
    Eigen::RowVectorXd Y = Eigen::RowVectorXd::Constant(N, 1.0);
    Eigen::RowVectorXd U = Eigen::RowVectorXd::Constant(N, 1.0);

    auto result = ctrlpp::batch_arx<2, 1>(Y, U);
    // Should not crash; system matrices should be finite
    REQUIRE(std::isfinite(result.system.A(0, 0)));
}

// ── MOESP hardening ────────────────────────────────────────────────────────────

TEST_CASE("MOESP with near-zero singular values", "[moesp][hardening][negative]")
{
    // Purely random noise data -- no underlying system
    constexpr int N = 100;
    Eigen::RowVectorXd Y(N);
    Eigen::RowVectorXd U(N);

    std::mt19937 gen(99);
    std::normal_distribution<double> dist(0.0, 0.001);

    for (int i = 0; i < N; ++i) {
        Y(i) = dist(gen);
        U(i) = dist(gen);
    }

    // Singular values should be near-zero, model won't be meaningful
    auto sv = ctrlpp::moesp_singular_values(Y, U);
    REQUIRE(sv.size() > 0);

    // Still identify a model -- should not crash
    auto result = ctrlpp::moesp<2>(Y, U);
    bool const cond_valid = std::isfinite(result.condition_number) || std::isinf(result.condition_number);
    REQUIRE(cond_valid);
}

TEST_CASE("MOESP degenerate data returns a zeroed system", "[moesp][hardening][negative]")
{
    // The documented failure contract (docs/api/sysid/moesp.md) is that a
    // failed identification reports condition_number = infinity and hands back
    // zeroed system matrices. discrete_state_space is an aggregate of Eigen
    // fixed-size matrices whose default constructor leaves coefficients
    // uninitialized, so omitting members yields indeterminate values unless the
    // zeroing is explicit -- this pins that it is.
    constexpr int N = 64;
    Eigen::RowVectorXd Y = Eigen::RowVectorXd::Zero(N);
    Eigen::RowVectorXd U = Eigen::RowVectorXd::Zero(N);

    auto result = ctrlpp::moesp<2>(Y, U);

    REQUIRE(std::isinf(result.condition_number));
    REQUIRE(result.system.A.isZero(0.0));
    REQUIRE(result.system.B.isZero(0.0));
    REQUIRE(result.system.C.isZero(0.0));
    REQUIRE(result.system.D.isZero(0.0));
}

TEST_CASE("MOESP with wrong model order", "[moesp][hardening][negative]")
{
    // True system is first order, identify with NX=4
    constexpr int N = 300;
    Eigen::RowVectorXd Y(N);
    Eigen::RowVectorXd U(N);

    std::mt19937 gen(123);
    std::uniform_real_distribution<double> input(-1.0, 1.0);
    std::normal_distribution<double> noise(0.0, 0.01);

    double x = 0.0;
    for (int t = 0; t < N; ++t) {
        U(t) = input(gen);
        Y(t) = x + noise(gen);
        x = 0.8 * x + 0.5 * U(t);
    }

    auto result = ctrlpp::moesp<4>(Y, U);
    REQUIRE(std::isfinite(result.system.A(0, 0)));
    // Over-specified model should still produce finite results
    REQUIRE(std::isfinite(result.metrics.nrmse));
}

TEST_CASE("MOESP identifies known state-space system", "[moesp][hardening][convergence]")
{
    // True system: x(t+1) = 0.9*x(t) + 0.5*u(t), y(t) = x(t)
    constexpr int N = 500;
    Eigen::RowVectorXd Y(N);
    Eigen::RowVectorXd U(N);

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> input(-1.0, 1.0);
    std::normal_distribution<double> noise(0.0, 0.01);

    double x = 0.0;
    for (int t = 0; t < N; ++t) {
        U(t) = input(gen);
        Y(t) = x + noise(gen);
        x = 0.9 * x + 0.5 * U(t);
    }

    auto result = ctrlpp::moesp<1>(Y, U);
    // NRMSE close to 0 indicates good fit (norm_error / norm_centered)
    REQUIRE(result.metrics.nrmse < 0.15);
}

TEST_CASE("Recursive ARX order 2 to_state_space superdiagonal",
          "[recursive_arx][hardening][coverage]")
{
    // NA=2 exercises the superdiagonal initialization loop in to_state_space()
    ctrlpp::recursive_arx<double, 2, 1> arx;

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> input(-1.0, 1.0);

    double y = 0.0, y_prev = 0.0;
    double u_prev = 0.0;
    for (int t = 0; t < 500; ++t) {
        double u = input(gen);
        double y_new = 0.6 * y + 0.2 * y_prev + 0.3 * u_prev;
        arx.update(y_new, u);
        y_prev = y;
        y = y_new;
        u_prev = u;
    }

    auto ss = arx.to_state_space();
    // 2nd-order system: A is 2x2 with superdiagonal entry
    REQUIRE(std::isfinite(ss.A(0, 1)));
    REQUIRE(ss.A.rows() == 2);
}
