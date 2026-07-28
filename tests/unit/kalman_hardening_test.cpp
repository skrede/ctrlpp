#include "hardening_helpers.h"
#include "ctrlpp/estimation/kalman.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

auto make_const_velocity_system()
{
    ctrlpp::discrete_state_space<double, 2, 1, 1> sys;
    double dt = 0.1;
    sys.A << 1.0, dt, 0.0, 1.0;
    sys.B << 0.5 * dt * dt, dt;
    sys.C << 1.0, 0.0;
    sys.D << 0.0;
    return sys;
}

}

TEST_CASE("Kalman rejects each non-finite configuration field by name", "[kalman][hardening][negative]")
{
    auto sys = make_const_velocity_system();

    // The failure this prevents: an infinite Q makes the first predict give
    // P = A P A' + Inf = Inf, the gain solve is posed against S = C P C' + R =
    // Inf and yields Inf/Inf = NaN, and the corrected state follows. The caller
    // would see a non-finite estimate from a filter it configured, with nothing
    // naming the field that was wrong.
    SECTION("process noise")
    {
        ctrlpp::kalman_config<double, 2, 1, 1> cfg{};
        cfg.Q = ctrlpp::test::inf_matrix<double, 2, 2>();
        const auto rejected = ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::non_finite_process_noise);
    }

    SECTION("measurement noise")
    {
        ctrlpp::kalman_config<double, 2, 1, 1> cfg{};
        cfg.R = ctrlpp::test::nan_matrix<double, 1, 1>();
        const auto rejected = ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::non_finite_measurement_noise);
    }

    SECTION("initial state")
    {
        ctrlpp::kalman_config<double, 2, 1, 1> cfg{};
        cfg.x0 = ctrlpp::test::nan_vector<double, 2>();
        const auto rejected = ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::non_finite_initial_state);
    }

    SECTION("initial covariance")
    {
        ctrlpp::kalman_config<double, 2, 1, 1> cfg{};
        cfg.P0 = ctrlpp::test::inf_matrix<double, 2, 2>();
        const auto rejected = ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::non_finite_initial_covariance);
    }
}

TEST_CASE("Kalman singular R (zero measurement noise)", "[kalman][hardening][negative]")
{
    auto sys = make_const_velocity_system();
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    Eigen::Matrix<double, 1, 1> R;
    R << 0.0;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Identity() * 10.0;

    auto kf = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}));

    // Should not crash -- singular R is a degenerate case
    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();
    kf.predict(u);

    Eigen::Matrix<double, 1, 1> z;
    z << 5.0;
    REQUIRE(kf.update(z).has_value());

    // State should be finite
    CHECK(std::isfinite(kf.state()[0]));
    CHECK(std::isfinite(kf.state()[1]));
}

TEST_CASE("Kalman NaN measurement is rejected without touching the estimate",
          "[kalman][hardening][negative]")
{
    auto sys = make_const_velocity_system();
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Identity();

    auto kf = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}));
    // Stepped only with the valid measurement, never with the poisoned one, so
    // it says what the filter would have carried had the bad sample never
    // arrived.
    auto reference = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}));

    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();
    kf.predict(u);
    reference.predict(u);

    // Snapshot immediately before the poisoned step.
    const Eigen::Vector2d x_before = kf.state();
    const Eigen::Matrix<double, 2, 2> P_before = kf.covariance();

    Eigen::Matrix<double, 1, 1> z_bad;
    z_bad << std::numeric_limits<double>::quiet_NaN();

    const auto rejected = kf.update(z_bad);

    REQUIRE_FALSE(rejected.has_value());
    REQUIRE(rejected.error() == ctrlpp::kalman_update_error::non_finite_measurement);

    // Exact comparison, not a tolerance: a rejected step performs no arithmetic
    // on the carried estimate at all, so bitwise equality is the contract and a
    // tolerance would admit a step that partially ran.
    CHECK(kf.state() == x_before);
    // The covariance was already measurement-independent before the guard
    // existed -- update_covariance(K) takes only the gain, never z -- so this
    // half of the invariant is structural. The state half is what the guard
    // adds.
    CHECK(kf.covariance() == P_before);
    // A rejection describes the sample, not the filter: nothing was mutated, so
    // the filter is not degraded and must not report that it is.
    CHECK(kf.health() == ctrlpp::kalman_health::ok);

    // The poison did not latch: the next valid step produces exactly what it
    // would have produced had the poisoned step never been attempted.
    Eigen::Matrix<double, 1, 1> z_good;
    z_good << 5.0;
    REQUIRE(kf.update(z_good).has_value());
    REQUIRE(reference.update(z_good).has_value());

    CHECK(kf.state() == reference.state());
    CHECK(kf.covariance() == reference.covariance());
}

TEST_CASE("Kalman zero Q", "[kalman][hardening][negative]")
{
    auto sys = make_const_velocity_system();
    auto Q = Eigen::Matrix<double, 2, 2>::Zero();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Identity();

    auto kf = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}));

    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();
    for(int k = 0; k < 100; ++k)
    {
        kf.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << 5.0;
        REQUIRE(kf.update(z).has_value());
    }

    CHECK(std::isfinite(kf.state()[0]));
    CHECK(std::isfinite(kf.state()[1]));
}

TEST_CASE("Kalman scalar analytical gain comparison", "[kalman][hardening][precision]")
{
    // Scalar system: A=1, B=0, C=1, D=0 (random walk)
    ctrlpp::discrete_state_space<double, 1, 1, 1> sys;
    sys.A << 1.0;
    sys.B << 0.0;
    sys.C << 1.0;
    sys.D << 0.0;

    double q_val = 0.1;
    double r_val = 1.0;
    double p0 = 10.0;

    Eigen::Matrix<double, 1, 1> Q, R, P0;
    Q << q_val;
    R << r_val;
    P0 << p0;
    Eigen::Matrix<double, 1, 1> x0;
    x0 << 0.0;

    auto kf = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 1, 1, 1>::create(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}));

    // Predict step: P_pred = A*P*A^T + Q = P + Q
    // After first predict: P_pred = p0 + q = 10.1
    // K = P_pred * H' / (H * P_pred * H' + R) = 10.1 / (10.1 + 1.0) = 10.1/11.1
    Eigen::Matrix<double, 1, 1> u_zero;
    u_zero << 0.0;
    kf.predict(u_zero);

    double P_pred = p0 + q_val;
    double expected_K = P_pred / (P_pred + r_val);

    Eigen::Matrix<double, 1, 1> z;
    z << 5.0;
    REQUIRE(kf.update(z).has_value());

    // After update: x = 0 + K*(5 - 0) = K*5
    REQUIRE_THAT(kf.state()[0], WithinAbs(expected_K * 5.0, 1e-12));
}

TEST_CASE("Kalman covariance stays positive definite over 1000 steps",
          "[kalman][hardening][stability]")
{
    auto sys = make_const_velocity_system();
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Identity() * 10.0;

    auto kf = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}));

    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();
    bool all_pd = true;

    for(int k = 0; k < 1000; ++k)
    {
        kf.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << 5.0 + 0.01 * k;
        REQUIRE(kf.update(z).has_value());

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> eigsolver(kf.covariance());
        for(int i = 0; i < 2; ++i)
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

TEST_CASE("Kalman state converges to truth within 200 steps", "[kalman][hardening][convergence]")
{
    auto sys = make_const_velocity_system();
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    Eigen::Matrix<double, 1, 1> R;
    R << 0.1;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Identity() * 10.0;

    auto kf = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}));

    double true_pos = 0.0;
    double true_vel = 1.0;
    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();
    double dt = 0.1;

    for(int k = 0; k < 200; ++k)
    {
        true_pos += true_vel * dt;
        kf.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << true_pos;
        REQUIRE(kf.update(z).has_value());
    }

    REQUIRE(std::abs(kf.state()[0] - true_pos) < 0.5);
}

TEST_CASE("Kalman ill-conditioned system matrix cond 1e10", "[kalman][hardening][robustness]")
{
    auto sys = ctrlpp::test::near_singular_system<double, 2, 1, 1>(1e-10);
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Identity();

    auto kf = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}));

    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();
    bool all_finite = true;

    for(int k = 0; k < 100; ++k)
    {
        kf.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << 1.0;
        REQUIRE(kf.update(z).has_value());

        if(!std::isfinite(kf.state()[0]) || !std::isfinite(kf.state()[1]))
        {
            all_finite = false;
            break;
        }
    }

    REQUIRE(all_finite);
}

TEST_CASE("Kalman is_steady_state with near-zero covariance",
          "[kalman][hardening][coverage]")
{
    constexpr std::size_t NX = 2, NU = 1, NY = 1;
    using KF = ctrlpp::kalman_filter<double, NX, NU, NY>;

    Eigen::Matrix2d A;
    A << 1.0, 0.01, 0.0, 1.0;
    Eigen::Vector2d B(0.0, 0.01);
    Eigen::RowVector2d C(1.0, 0.0);
    Eigen::Matrix<double, 1, 1> D = Eigen::Matrix<double, 1, 1>::Zero();
    ctrlpp::discrete_state_space<double, NX, NU, NY> sys{A, B, C, D};

    ctrlpp::kalman_config<double, NX, NU, NY> cfg{};
    // Near-zero initial covariance triggers P.norm() < epsilon branch
    cfg.P0 = Eigen::Matrix2d::Zero();
    cfg.Q = Eigen::Matrix2d::Identity() * 1e-300;
    cfg.R = Eigen::Matrix<double, 1, 1>::Identity();

    auto kf = ctrlpp::test::constructed(KF::create(sys, cfg));

    // With zero P, should immediately report steady state
    CHECK(kf.is_steady_state());
}
