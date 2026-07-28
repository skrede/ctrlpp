#include "hardening_helpers.h"
#include "ctrlpp/estimation/ukf.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;

namespace {

struct ukf_linear_dynamics
{
    double dt = 0.1;

    auto operator()(const ctrlpp::Vector<double, 2>& x,
                    const ctrlpp::Vector<double, 1>& u) const -> ctrlpp::Vector<double, 2>
    {
        ctrlpp::Vector<double, 2> xn;
        xn(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        xn(1) = x(1) + dt * u(0);
        return xn;
    }
};

struct ukf_position_measurement
{
    auto operator()(const ctrlpp::Vector<double, 2>& x) const -> ctrlpp::Vector<double, 1>
    {
        ctrlpp::Vector<double, 1> z;
        z(0) = x(0);
        return z;
    }
};

auto make_ukf()
{
    ctrlpp::ukf_config<double, 2, 1, 1> cfg{};
    cfg.Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    cfg.R = Eigen::Matrix<double, 1, 1>::Identity();
    cfg.x0 = Eigen::Vector2d::Zero();
    cfg.P0 = Eigen::Matrix<double, 2, 2>::Identity() * 10.0;
    return ctrlpp::test::constructed(ctrlpp::ukf<double, 2, 1, 1, ukf_linear_dynamics, ukf_position_measurement>::create(ukf_linear_dynamics{}, ukf_position_measurement{}, cfg));
}

}

TEST_CASE("UKF rejects each non-finite configuration field by name", "[ukf][hardening][negative]")
{
    using filter_t = ctrlpp::ukf<double, 2, 1, 1, ukf_linear_dynamics, ukf_position_measurement>;

    // The failure this prevents: an infinite Q is added to the sigma-point
    // covariance every predict, so the spread the next generation factors is
    // infinite, the gain solve is posed against an infinite innovation
    // covariance, and the estimate is non-finite at the first step.
    SECTION("process noise")
    {
        ctrlpp::ukf_config<double, 2, 1, 1> cfg{};
        cfg.Q = ctrlpp::test::inf_matrix<double, 2, 2>();
        const auto rejected = filter_t::create(ukf_linear_dynamics{}, ukf_position_measurement{}, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::non_finite_process_noise);
    }

    SECTION("measurement noise")
    {
        ctrlpp::ukf_config<double, 2, 1, 1> cfg{};
        cfg.R = ctrlpp::test::nan_matrix<double, 1, 1>();
        const auto rejected = filter_t::create(ukf_linear_dynamics{}, ukf_position_measurement{}, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::non_finite_measurement_noise);
    }

    SECTION("initial state")
    {
        ctrlpp::ukf_config<double, 2, 1, 1> cfg{};
        cfg.x0 = ctrlpp::test::nan_vector<double, 2>();
        const auto rejected = filter_t::create(ukf_linear_dynamics{}, ukf_position_measurement{}, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::non_finite_initial_state);
    }

    SECTION("initial covariance")
    {
        ctrlpp::ukf_config<double, 2, 1, 1> cfg{};
        cfg.P0 = ctrlpp::test::inf_matrix<double, 2, 2>();
        const auto rejected = filter_t::create(ukf_linear_dynamics{}, ukf_position_measurement{}, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::non_finite_initial_covariance);
    }

    // The options-taking overload reports the strategy's own rejection FIRST: a
    // strategy that cannot be built leaves nothing for the configuration to be a
    // configuration of. Both faults are present here and the strategy wins.
    SECTION("a rejected strategy is reported ahead of a rejected configuration")
    {
        ctrlpp::ukf_config<double, 2, 1, 1> cfg{};
        cfg.Q = ctrlpp::test::inf_matrix<double, 2, 2>();
        const auto rejected = filter_t::create(ukf_linear_dynamics{}, ukf_position_measurement{}, cfg,
                                               ctrlpp::merwe_options<double>{.alpha = 0.0, .beta = 2.0, .kappa = 0.0});
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::non_positive_sigma_spread);
    }
}

TEST_CASE("UKF NaN measurement is rejected without touching the estimate",
          "[ukf][hardening][negative]")
{
    auto filter = make_ukf();
    // Stepped only with the valid measurement, never with the poisoned one, so
    // it says what the filter would have carried had the bad sample never
    // arrived.
    auto reference = make_ukf();

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;
    filter.predict(u);
    reference.predict(u);

    // Snapshot immediately before the poisoned step.
    const Eigen::Vector2d x_before = filter.state();
    const Eigen::Matrix<double, 2, 2> P_before = filter.covariance();

    Eigen::Matrix<double, 1, 1> z_bad;
    z_bad << std::numeric_limits<double>::quiet_NaN();

    const auto rejected = filter.update(z_bad);

    REQUIRE_FALSE(rejected.has_value());
    REQUIRE(rejected.error() == ctrlpp::ukf_update_error::non_finite_measurement);

    // Exact comparison, not a tolerance: a rejected step performs no arithmetic
    // on the carried estimate at all, so bitwise equality is the contract and a
    // tolerance would admit a step that partially ran.
    CHECK(filter.state() == x_before);
    // The covariance was already measurement-independent before the guard
    // existed -- the posterior reduction P - K*S*K^T is built from the
    // predict-stage sigma points and the gain, never from z -- so this half of
    // the invariant is structural. The state half is what the guard adds.
    CHECK(filter.covariance() == P_before);
    // A rejection describes the sample, not the filter: nothing was mutated, so
    // the filter is not degraded and must not report that it is.
    CHECK(filter.health() == ctrlpp::ukf_health::ok);

    // The poison did not latch: the next valid step produces exactly what it
    // would have produced had the poisoned step never been attempted.
    Eigen::Matrix<double, 1, 1> z_good;
    z_good << 1.0;
    REQUIRE(filter.update(z_good).has_value());
    REQUIRE(reference.update(z_good).has_value());

    CHECK(filter.state() == reference.state());
    CHECK(filter.covariance() == reference.covariance());
}

TEST_CASE("UKF for linear system matches Kalman output", "[ukf][hardening][precision]")
{
    auto filter = make_ukf();

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;

    double true_pos = 0.0;
    double true_vel = 1.0;

    for(int k = 0; k < 50; ++k)
    {
        true_pos += 0.1 * true_vel;
        filter.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << true_pos;
        REQUIRE(filter.update(z).has_value());
    }

    REQUIRE(std::abs(filter.state()[0] - true_pos) < 1.0);
}

TEST_CASE("UKF covariance stays PD over 1000 steps", "[ukf][hardening][stability]")
{
    auto filter = make_ukf();

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;
    bool all_pd = true;

    for(int k = 0; k < 1000; ++k)
    {
        filter.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << 5.0 + 0.01 * k;
        REQUIRE(filter.update(z).has_value());

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> eigsolver(filter.covariance());
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

TEST_CASE("UKF tracks nonlinear system (quadratic dynamics)", "[ukf][hardening][convergence]")
{
    struct quadratic_dynamics
    {
        auto operator()(const ctrlpp::Vector<double, 2>& x,
                        const ctrlpp::Vector<double, 1>& /*u*/) const -> ctrlpp::Vector<double, 2>
        {
            ctrlpp::Vector<double, 2> xn;
            xn(0) = 0.95 * x(0) + 0.05 * x(1);
            xn(1) = 0.95 * x(1);
            return xn;
        }
    };

    ctrlpp::ukf_config<double, 2, 1, 1> cfg{};
    cfg.Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    cfg.R = Eigen::Matrix<double, 1, 1>::Identity() * 0.1;
    cfg.x0 = Eigen::Vector2d::Zero();
    cfg.P0 = Eigen::Matrix<double, 2, 2>::Identity() * 10.0;

    auto filter = ctrlpp::test::constructed(ctrlpp::ukf<double, 2, 1, 1, quadratic_dynamics, ukf_position_measurement>::create(quadratic_dynamics{}, ukf_position_measurement{}, cfg));

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;

    // True state converges to zero -- filter should track
    ctrlpp::Vector<double, 2> x_true;
    x_true << 5.0, 1.0;

    for(int k = 0; k < 200; ++k)
    {
        x_true(0) = 0.95 * x_true(0) + 0.05 * x_true(1);
        x_true(1) = 0.95 * x_true(1);

        filter.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << x_true(0);
        REQUIRE(filter.update(z).has_value());
    }

    REQUIRE(std::abs(filter.state()[0] - x_true(0)) < 0.5);
}

TEST_CASE("UKF near-singular P0", "[ukf][hardening][robustness]")
{
    ctrlpp::ukf_config<double, 2, 1, 1> cfg{};
    cfg.Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    cfg.R = Eigen::Matrix<double, 1, 1>::Identity();
    cfg.x0 = Eigen::Vector2d::Zero();
    cfg.P0 = ctrlpp::test::ill_conditioned_2x2<double>(1e10);

    // Deliberately near-singular and entirely well-posed: every entry is finite,
    // so the configuration validation accepts it, and whether the Cholesky-based
    // sigma-point generation survives it is the question the case exists to ask.
    // A validation that rejected a small-but-positive covariance entry would
    // delete that question.
    auto filter = ctrlpp::test::constructed(ctrlpp::ukf<double, 2, 1, 1, ukf_linear_dynamics, ukf_position_measurement>::create(ukf_linear_dynamics{}, ukf_position_measurement{}, cfg));

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;

    bool all_finite = true;
    for(int k = 0; k < 100; ++k)
    {
        filter.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << 1.0;
        REQUIRE(filter.update(z).has_value());

        if(!std::isfinite(filter.state()[0]) || !std::isfinite(filter.state()[1]))
        {
            all_finite = false;
            break;
        }
    }

    REQUIRE(all_finite);
}
