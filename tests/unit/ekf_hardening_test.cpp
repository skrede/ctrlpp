#include "hardening_helpers.h"
#include "ctrlpp/estimation/ekf.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

struct linear_dynamics
{
    auto operator()(const ctrlpp::Vector<double, 2>& x,
                    const ctrlpp::Vector<double, 1>& u) const -> ctrlpp::Vector<double, 2>
    {
        ctrlpp::Vector<double, 2> xn;
        xn[0] = x[0] + 0.1 * x[1];
        xn[1] = x[1] + 0.1 * u[0];
        return xn;
    }
};

struct linear_measurement
{
    auto operator()(const ctrlpp::Vector<double, 2>& x) const -> ctrlpp::Vector<double, 1>
    {
        ctrlpp::Vector<double, 1> z;
        z << x[0];
        return z;
    }
};

using ekf_t = ctrlpp::ekf<double, 2, 1, 1, linear_dynamics, linear_measurement>;

// create() is the only construction path and it is fallible, so every
// valid-input site goes through it and asserts success here. Cases that mean to
// observe a rejected configuration assert on the result directly instead.
auto build_ekf(const ctrlpp::ekf_config<double, 2, 1, 1>& cfg) -> ekf_t
{
    return ctrlpp::test::constructed(ekf_t::create(linear_dynamics{}, linear_measurement{}, cfg));
}

auto make_ekf()
{
    ctrlpp::ekf_config<double, 2, 1, 1> cfg{};
    cfg.Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    cfg.R = Eigen::Matrix<double, 1, 1>::Identity();
    cfg.x0 = Eigen::Vector2d::Zero();
    cfg.P0 = Eigen::Matrix<double, 2, 2>::Identity() * 10.0;

    return build_ekf(cfg);
}

}

TEST_CASE("EKF NaN measurement is rejected without touching the estimate",
          "[ekf][hardening][negative]")
{
    auto filter = make_ekf();
    // Stepped only with the valid measurement, never with the poisoned one, so
    // it says what the filter would have carried had the bad sample never
    // arrived.
    auto reference = make_ekf();

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
    REQUIRE(rejected.error() == ctrlpp::ekf_update_error::non_finite_measurement);

    // Exact comparison, not a tolerance: a rejected step performs no arithmetic
    // on the carried estimate at all, so bitwise equality is the contract and a
    // tolerance would admit a step that partially ran.
    CHECK(filter.state() == x_before);
    // The covariance was already measurement-independent before the guard
    // existed -- update_covariance(K, H) takes the gain and the measurement
    // Jacobian, never z -- so this half of the invariant is structural. The
    // state half is what the guard adds.
    CHECK(filter.covariance() == P_before);
    // A rejection describes the sample, not the filter: nothing was mutated, so
    // the filter is not degraded and must not report that it is.
    CHECK(filter.health() == ctrlpp::ekf_health::ok);

    // The poison did not latch: the next valid step produces exactly what it
    // would have produced had the poisoned step never been attempted.
    Eigen::Matrix<double, 1, 1> z_good;
    z_good << 1.0;
    REQUIRE(filter.update(z_good).has_value());
    REQUIRE(reference.update(z_good).has_value());

    CHECK(filter.state() == reference.state());
    CHECK(filter.covariance() == reference.covariance());
}

TEST_CASE("EKF infinite process noise is rejected at construction", "[ekf][hardening][negative]")
{
    ctrlpp::ekf_config<double, 2, 1, 1> cfg{};
    cfg.Q = Eigen::Matrix<double, 2, 2>::Identity() * std::numeric_limits<double>::infinity();
    cfg.R = Eigen::Matrix<double, 1, 1>::Identity();
    cfg.x0 = Eigen::Vector2d::Zero();
    cfg.P0 = Eigen::Matrix<double, 2, 2>::Identity();

    // Had the configuration been accepted, the fault would have surfaced as far
    // as possible from where it was made. Q is added to the propagated
    // covariance, so the first predict gives P = F P F' + Inf = Inf; the gain
    // solve is then posed against S = H P H' + R = Inf and yields Inf/Inf, which
    // is NaN in every entry; the corrected state x + K y is NaN with it, and the
    // Joseph-form covariance follows. The caller would see a non-finite estimate
    // from a filter it configured and would have no way to tell which field was
    // wrong. Rejecting here names the field instead.
    const auto rejected = ekf_t::create(linear_dynamics{}, linear_measurement{}, cfg);

    REQUIRE_FALSE(rejected.has_value());
    REQUIRE(rejected.error() == ctrlpp::filter_error::non_finite_process_noise);
}

TEST_CASE("EKF rejects each non-finite configuration field by name", "[ekf][hardening][negative]")
{
    const auto inf2 = ctrlpp::test::inf_matrix<double, 2, 2>();

    SECTION("measurement noise")
    {
        ctrlpp::ekf_config<double, 2, 1, 1> cfg{};
        cfg.R = ctrlpp::test::nan_matrix<double, 1, 1>();
        const auto rejected = ekf_t::create(linear_dynamics{}, linear_measurement{}, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::non_finite_measurement_noise);
    }

    SECTION("initial state")
    {
        ctrlpp::ekf_config<double, 2, 1, 1> cfg{};
        cfg.x0 = ctrlpp::test::nan_vector<double, 2>();
        const auto rejected = ekf_t::create(linear_dynamics{}, linear_measurement{}, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::non_finite_initial_state);
    }

    SECTION("initial covariance")
    {
        ctrlpp::ekf_config<double, 2, 1, 1> cfg{};
        cfg.P0 = inf2;
        const auto rejected = ekf_t::create(linear_dynamics{}, linear_measurement{}, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::non_finite_initial_covariance);
    }
}

TEST_CASE("EKF for linear system matches Kalman gain within 1%", "[ekf][hardening][precision]")
{
    // For a linear system, the EKF should produce the same result as the Kalman filter
    auto filter = make_ekf();

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

    // After 50 steps, position estimate should be close to truth
    REQUIRE(std::abs(filter.state()[0] - true_pos) < 1.0);
}

TEST_CASE("EKF covariance stays PD over 1000 steps", "[ekf][hardening][stability]")
{
    auto filter = make_ekf();

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

TEST_CASE("EKF state estimate converges to truth", "[ekf][hardening][convergence]")
{
    auto filter = make_ekf();

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;

    double true_pos = 0.0;
    double true_vel = 1.0;

    for(int k = 0; k < 200; ++k)
    {
        true_pos += 0.1 * true_vel;
        filter.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << true_pos;
        REQUIRE(filter.update(z).has_value());
    }

    REQUIRE(std::abs(filter.state()[0] - true_pos) < 0.5);
}

TEST_CASE("EKF ill-conditioned system cond 1e10", "[ekf][hardening][robustness]")
{
    ctrlpp::ekf_config<double, 2, 1, 1> cfg{};
    cfg.Q = ctrlpp::test::ill_conditioned_2x2<double>(1e10);
    cfg.R = Eigen::Matrix<double, 1, 1>::Identity();
    cfg.x0 = Eigen::Vector2d::Zero();
    cfg.P0 = Eigen::Matrix<double, 2, 2>::Identity();

    // Deliberately ill-conditioned and entirely well-posed: every entry is
    // finite, so the configuration validation accepts it. Conditioning is a
    // numerical-behavior question and finiteness is the domain condition; a
    // validation that rejected this would refuse a problem the filter solves.
    auto filter = build_ekf(cfg);

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
