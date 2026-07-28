// What the oracles in this file decide.
//
// Every dynamics model in this file is AFFINE, and the unscented transform is
// exact for an affine map: the propagated sigma set recombines to exactly the
// image of the mean and exactly the congruence of the covariance. So the
// reference for every tracking case here is an independently constructed
// ctrlpp::kalman_filter on the same data, and the only thing separating the two
// filters is arithmetic.
//
// That arithmetic is the interesting part, and the budgets say so. The scaled
// sigma-point weights with the default spread are of magnitude
// |lambda| / (n + lambda), which is about a million: the recombination sums
// terms a million times larger than their own total, so the transform's answer
// carries a million epsilons of cancellation. Every budget below is that
// amplification -- computed in the test from the strategy's own parameters, not
// written down as a number -- times a counted operation total times the operand
// scale. The consequence is worth stating plainly: this filter agrees with the
// linear one to about a part in 1e9, not to a part in 1e16, and the reason is
// the spread parameter rather than anything about the model.
//
// What they deliberately do not decide. Nothing here asserts nonlinear tracking
// accuracy. The case that used to claim it named its dynamics quadratic and
// supplied affine ones, so it never tested nonlinearity at all; it is renamed
// below rather than replaced, because a genuinely nonlinear model has no
// closed-form reference and substituting a second implementation for an oracle
// would restore exactly the defect this file is being cleared of.

#include "hardening_helpers.h"

#include "ctrlpp/estimation/ukf.h"
#include "ctrlpp/estimation/kalman.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;

namespace {

constexpr double ukf_eps = std::numeric_limits<double>::epsilon();

// Rounded operations at the sigma-point weights' own magnitude, along one
// filter step. Five weighted recombinations occur -- the predicted mean, the
// predicted covariance, the predicted measurement, the innovation covariance
// and the cross covariance -- and each sums 2n+1 = 5 terms, so 25 operations
// carry the cancellation. The gain solve and the correction add ordinary
// roundings on top of that and are dominated by it.
constexpr int ukf_transform_ops = 25;

// The symmetric eigensolver's backward perturbation, bounded by the number of
// entries of a two-by-two times epsilon times the matrix norm (Weyl).
constexpr int ukf_eig_ops = 4;

/// @brief Cancellation amplification of the scaled unscented transform, derived
/// from the strategy's own parameters rather than measured after the fact.
///
/// The scaling term is lambda = alpha^2 (n + kappa) - n and the weight
/// denominator is n + lambda = alpha^2 (n + kappa), so the centre weight is
/// lambda / (n + lambda) and the outer ones are half its reciprocal denominator.
/// The recombination is a sum of terms of that magnitude whose total is of unit
/// magnitude, so the largest weight IS the factor by which a relative rounding
/// in a term becomes a relative error in the result.
auto unscented_amplification(const ctrlpp::merwe_options<double>& opts, double n) -> double
{
    const double lambda = opts.alpha * opts.alpha * (n + opts.kappa) - n;
    return std::abs(lambda) / (n + lambda);
}

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

/// @brief The linear system `ukf_linear_dynamics` and the position measurement
/// together ARE, written out so an independent linear filter can be driven with
/// the same data.
auto equivalent_system() -> ctrlpp::discrete_state_space<double, 2, 1, 1>
{
    ctrlpp::discrete_state_space<double, 2, 1, 1> sys;
    const double dt = ukf_linear_dynamics{}.dt;
    sys.A << 1.0, dt, 0.0, 1.0;
    sys.B << 0.5 * dt * dt, dt;
    sys.C << 1.0, 0.0;
    sys.D << 0.0;
    return sys;
}

/// @brief Spectral norm of the closed-loop error map the reference filter's own
/// reported covariance implies. It exceeds one during the transient, so a bound
/// formed as a plain multiple of the step count would not be a bound.
auto closed_loop_norm(const ctrlpp::discrete_state_space<double, 2, 1, 1>& sys,
                      const Eigen::Matrix2d& P_pred, double R) -> double
{
    const double S = (sys.C * P_pred * sys.C.transpose())(0, 0) + R;
    const Eigen::Vector2d K = (P_pred * sys.C.transpose()).eval() / S;
    return ((Eigen::Matrix2d::Identity() - K * sys.C) * sys.A).norm();
}

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

TEST_CASE("UKF for a linear system IS the Kalman filter", "[ukf][hardening][precision]")
{
    // The name this case used to carry promised agreement with the Kalman filter
    // and then constructed none: it asserted the position estimate was within 1.0
    // of the truth, which is one measurement standard deviation written as a bare
    // literal and true of a filter that ignored its measurements.
    const auto sys = equivalent_system();
    const double amplitude = unscented_amplification(ctrlpp::merwe_options<double>{}, 2.0);

    auto filter = make_ukf();
    auto reference = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(
        sys, {.Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01,
              .R = Eigen::Matrix<double, 1, 1>::Identity(),
              .x0 = Eigen::Vector2d::Zero(),
              .P0 = Eigen::Matrix<double, 2, 2>::Identity() * 10.0}));

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;

    double true_pos = 0.0;
    double divergence_budget = 0.0;
    double worst_divergence = 0.0;
    double worst_scale = 1.0;

    for(int k = 0; k < 50; ++k)
    {
        CAPTURE(k);
        true_pos += 0.1;
        filter.predict(u);
        reference.predict(u);

        const double propagation = closed_loop_norm(sys, reference.covariance(), 1.0);

        Eigen::Matrix<double, 1, 1> z;
        z << true_pos;
        REQUIRE(filter.update(z).has_value());
        REQUIRE(reference.update(z).has_value());

        const double scale = std::max({filter.state().cwiseAbs().maxCoeff(),
                                       reference.state().cwiseAbs().maxCoeff(), 1.0});
        worst_scale = std::max(worst_scale, scale);
        divergence_budget =
            propagation * divergence_budget + ukf_transform_ops * amplitude * ukf_eps * scale;

        CAPTURE(filter.state()[0], reference.state()[0], divergence_budget);
        REQUIRE_THAT(filter.state()[0], WithinAbs(reference.state()[0], divergence_budget));
        REQUIRE_THAT(filter.state()[1], WithinAbs(reference.state()[1], divergence_budget));

        worst_divergence =
            std::max(worst_divergence, (filter.state() - reference.state()).cwiseAbs().maxCoeff());
    }

    // The amplification is the budget's only factor that is not a counted
    // operation, so the case pins it from the other side too: the realized
    // disagreement is well ABOVE what the same operation count would allow
    // without it. That makes the amplification load-bearing rather than
    // decorative -- were the default spread changed to one, the transform would
    // agree with the linear filter to plain arithmetic rounding, this assertion
    // would fail, and every budget in this file would need revisiting. Which is
    // the correct outcome, not a maintenance hazard.
    CAPTURE(worst_divergence, amplitude);
    REQUIRE(worst_divergence > ukf_transform_ops * ukf_eps * worst_scale);
}

TEST_CASE("UKF covariance stays PD over 1000 steps", "[ukf][hardening][stability]")
{
    auto filter = make_ukf();

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;

    for(int k = 0; k < 1000; ++k)
    {
        CAPTURE(k);
        filter.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << 5.0 + 0.01 * k;
        REQUIRE(filter.update(z).has_value());

        // Symmetry is the other half of the covariance contract and was checked
        // nowhere. An equality, not a tolerance: the filter symmetrizes both the
        // predicted and the corrected covariance.
        REQUIRE(filter.covariance()(0, 1) == filter.covariance()(1, 0));

        // The floor is the eigensolver's own backward error at the scale of the
        // matrix handed to it. The threshold it replaces, -1e-10, is about
        // 450000 epsilons of slack BELOW zero.
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> eigsolver(filter.covariance());
        REQUIRE(eigsolver.eigenvalues().minCoeff()
                >= -ukf_eig_ops * ukf_eps * filter.covariance().norm());
    }

    // Every covariance factored over a thousand steps was positive definite on
    // its own terms: no repair was needed, and the filter says so. A repaired
    // covariance is a degraded regime the health status exists to surface, so
    // asserting it never happened is a stronger statement than the eigenvalue
    // floor alone -- the floor tolerates a matrix that was repaired INTO
    // definiteness, and this does not.
    REQUIRE(filter.health() == ctrlpp::ukf_health::ok);
}

TEST_CASE("UKF tracks a decaying coupled system exactly as the linear filter does",
          "[ukf][hardening][convergence]")
{
    // This case used to be named for quadratic dynamics and supplied these, which
    // are AFFINE in both coordinates: every term is a constant times a state
    // component and there is no product of components anywhere. The case
    // therefore never exercised nonlinearity, and its bound of 0.5 was fitted
    // against a realized error of order 1e-9. Naming what the dynamics are and
    // asserting the equality they imply is the honest repair; the alternative,
    // making them genuinely nonlinear, would leave the case with no reference to
    // be checked against.
    struct decaying_coupled_dynamics
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

    auto filter = ctrlpp::test::constructed(ctrlpp::ukf<double, 2, 1, 1, decaying_coupled_dynamics, ukf_position_measurement>::create(decaying_coupled_dynamics{}, ukf_position_measurement{}, cfg));

    ctrlpp::discrete_state_space<double, 2, 1, 1> sys;
    sys.A << 0.95, 0.05, 0.0, 0.95;
    sys.B << 0.0, 0.0;
    sys.C << 1.0, 0.0;
    sys.D << 0.0;
    auto reference = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(
        sys, {.Q = cfg.Q, .R = cfg.R, .x0 = cfg.x0, .P0 = cfg.P0}));

    const double amplitude = unscented_amplification(ctrlpp::merwe_options<double>{}, 2.0);

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;

    ctrlpp::Vector<double, 2> x_true;
    x_true << 5.0, 1.0;
    double divergence_budget = 0.0;

    for(int k = 0; k < 200; ++k)
    {
        CAPTURE(k);
        x_true = (sys.A * x_true).eval();

        filter.predict(u);
        reference.predict(u);

        const double propagation = closed_loop_norm(sys, reference.covariance(), cfg.R(0, 0));

        Eigen::Matrix<double, 1, 1> z;
        z << x_true(0);
        REQUIRE(filter.update(z).has_value());
        REQUIRE(reference.update(z).has_value());

        const double scale = std::max({x_true.cwiseAbs().maxCoeff(),
                                       filter.state().cwiseAbs().maxCoeff(),
                                       reference.state().cwiseAbs().maxCoeff(), 1.0});
        divergence_budget =
            propagation * divergence_budget + ukf_transform_ops * amplitude * ukf_eps * scale;

        CAPTURE(filter.state()[0], reference.state()[0], x_true(0));
        REQUIRE_THAT(filter.state()[0], WithinAbs(reference.state()[0], divergence_budget));
        REQUIRE_THAT(filter.state()[1], WithinAbs(reference.state()[1], divergence_budget));
    }
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

    // The reference is the linear filter on the identical near-singular data.
    // Whether the Cholesky-based sigma-point generation survives is exactly what
    // finiteness could not see: a factorization that had silently degraded would
    // still produce finite numbers and would immediately part company with the
    // reference.
    const auto sys = equivalent_system();
    auto reference = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(
        sys, {.Q = cfg.Q, .R = cfg.R, .x0 = cfg.x0, .P0 = cfg.P0}));

    const double amplitude = unscented_amplification(ctrlpp::merwe_options<double>{}, 2.0);

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;
    double divergence_budget = 0.0;

    for(int k = 0; k < 100; ++k)
    {
        CAPTURE(k);
        filter.predict(u);
        reference.predict(u);

        const double propagation = closed_loop_norm(sys, reference.covariance(), cfg.R(0, 0));

        Eigen::Matrix<double, 1, 1> z;
        z << 1.0;
        REQUIRE(filter.update(z).has_value());
        REQUIRE(reference.update(z).has_value());

        const double scale = std::max({filter.state().cwiseAbs().maxCoeff(),
                                       reference.state().cwiseAbs().maxCoeff(), 1.0});
        divergence_budget =
            propagation * divergence_budget + ukf_transform_ops * amplitude * ukf_eps * scale;

        CAPTURE(filter.state()[0], reference.state()[0]);
        REQUIRE_THAT(filter.state()[0], WithinAbs(reference.state()[0], divergence_budget));
        REQUIRE_THAT(filter.state()[1], WithinAbs(reference.state()[1], divergence_budget));

        REQUIRE(filter.covariance()(0, 1) == filter.covariance()(1, 0));
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> eigsolver(filter.covariance());
        REQUIRE(eigsolver.eigenvalues().minCoeff()
                >= -ukf_eig_ops * ukf_eps * filter.covariance().norm());
    }

    // The near-singular initial covariance never forced a repair: the factorization
    // handled it on its own terms, and the filter says so. That is the case's
    // actual question, and it now has an answer.
    REQUIRE(filter.health() == ctrlpp::ukf_health::ok);

    // The measurement is constant, so the observed coordinate must have settled on
    // it to within the uncertainty the filter itself reports. A multiplier of one
    // is the consistency contract: an error above the reported standard deviation
    // means the filter claims a confidence its estimate does not earn.
    CAPTURE(filter.state()[0], filter.covariance()(0, 0));
    REQUIRE(std::abs(filter.state()[0] - 1.0) <= std::sqrt(filter.covariance()(0, 0)));
}
