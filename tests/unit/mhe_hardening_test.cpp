// What the oracles in this file decide.
//
// The estimator reports on TWO channels and this file asserts both. A refusal
// -- "I could not do what you asked" -- is the return value of `update`. A
// disposition on a step that succeeded -- which of the two estimators produced
// the estimate, and how the solve went -- is the diagnostics aggregate. No case
// in this file read either before:
//
//  * A poisoned measurement is REFUSED by the embedded filter, and the estimator
//    forwards that verdict verbatim on the failure channel, naming the operand.
//    It invalidates the fixed-step window but keeps the prediction for the
//    input the plant already received. The diagnostics aggregate is unchanged
//    across the refusal because it describes the last successful update. The
//    estimator then uses the embedded filter until a coherent horizon refills.
//  * During the fill-up the estimator delegates to its embedded filter, and it
//    says so: the fallback flag is true for exactly the first N updates and
//    false from the one that first fills the window. That is the documented
//    warm-up contract and nothing tested it.
//  * Grossly inconsistent measurements are NOT a failure. The window solve
//    still reports an optimal status with the fallback cleared, and the
//    estimate it produces is exactly ODD in the record: negating every
//    measurement negates the estimate BITWISE, which is a real property of a
//    quadratic cost with no offset and which no fitted bound could state.
//  * Agreement with a linear filter is asserted by CONSTRUCTING that filter --
//    the same dynamics, the same weights, the same initial condition -- and
//    comparing state for state at every step. The identity is inexact and the
//    inexactness is pinned rather than hidden: the arrival cost is the embedded
//    filter's covariance at the window head rather than the exact conditional
//    one, so the two estimators disagree at the level the solver's own
//    tolerance permits.
//  * Convergence is asserted against the estimator's OWN reported covariance,
//    at three standard errors, which is the convention this phase's estimation
//    oracles already use -- not against a count of steps whose error happened
//    to fall, which is a statistic about one seed rather than a contract.
//  * An ill-conditioned process weight is asserted through what it MEANS: the
//    two process variances differ by ten decades, so the penalties on moving
//    each coordinate differ by ten decades, so the velocity may travel across
//    the window only that fraction of what the position travels. And the solver
//    reports `solved_inaccurate` rather than `optimal`, which is the honest
//    disposition for this data and which no case observed.
//
// What they deliberately do not decide. Nothing here reconstructs the window
// least-squares solution in closed form. Doing so would be a second
// implementation of the estimator under test, which is the defect class this
// phase exists to remove wearing different clothes; the odd symmetry, the
// filter comparison and the frozen velocity state real properties without it.

#include "hardening_helpers.h"

#include "ctrlpp/mhe.h"
#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/estimation/kalman.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>

#include <cmath>
#include <cstddef>
#include <limits>
#include <random>

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

// The conventional bound on an estimate against its own reported standard
// deviation, stated once and used unchanged wherever a statistical claim is
// made here. It is not an operation count.
constexpr double reported_standard_errors = 3.0;

}

// ── MHE hardening: negative ────────────────────────────────────────────────────

TEST_CASE("MHE with NaN in measurement noise", "[mhe][hardening][negative]")
{
    ctrlpp::mhe_config<double, NX, NU, NY, N> cfg;
    cfg.Q = ctrlpp::Matrix<double, NX, NX>::Identity() * 0.01;
    cfg.R = ctrlpp::Matrix<double, NY, NY>::Identity() * 0.1;

    auto estimator = ctrlpp::test::constructed(MheType::create(linear_dynamics{}, position_measurement{}, cfg));

    ctrlpp::Vector<double, NU> u = ctrlpp::Vector<double, NU>::Zero();

    // Warm up with valid measurements, and assert the documented warm-up
    // contract while doing it: the estimator delegates to its embedded filter
    // until the window is full, and the aggregate says which of the two
    // produced the estimate.
    for(std::size_t i = 0; i < N + 1; ++i)
    {
        estimator.predict(u);
        ctrlpp::Vector<double, NY> z;
        z << 0.1 * static_cast<double>(i);
        REQUIRE(estimator.update(z).has_value());

        REQUIRE(estimator.diagnostics().status == ctrlpp::solve_status::optimal);
        REQUIRE(estimator.diagnostics().used_ekf_fallback == (i < N));
    }
    REQUIRE(estimator.is_initialized());

    const ctrlpp::solve_status status_before = estimator.diagnostics().status;
    const bool fallback_before = estimator.diagnostics().used_ekf_fallback;
    const double cost_before = estimator.diagnostics().cost;

    // Inject NaN measurement. The covariance is snapshotted AFTER the
    // propagation, because the reported covariance is the embedded filter's and
    // the propagation legitimately advances it; the claim is about the update.
    estimator.predict(u);
    const ctrlpp::Matrix<double, NX, NX> covariance_before = estimator.covariance();
    ctrlpp::Vector<double, NY> z_nan;
    z_nan << std::numeric_limits<double>::quiet_NaN();

    const auto refused = estimator.update(z_nan);

    // The failure channel, which is what a refusal belongs on: the caller is
    // told that the step did not happen and told which operand caused it,
    // forwarded verbatim from the embedded filter rather than renamed.
    REQUIRE_FALSE(refused.has_value());
    CHECK(refused.error() == ctrlpp::ekf_update_error::non_finite_measurement);

    // The input was already applied and its prediction remains current, but the
    // fixed-step window cannot express the missing measurement and is reset.
    CHECK(estimator.state().allFinite());
    CHECK(estimator.covariance() == covariance_before);
    CHECK_FALSE(estimator.is_initialized());

    // And the disposition channel says nothing about it, which is the point of
    // having two channels: it still describes the last step that SUCCEEDED, the
    // same step `state()` describes. A refusal is not a disposition on a
    // result, because there is no result.
    CHECK(estimator.diagnostics().status == status_before);
    CHECK(estimator.diagnostics().used_ekf_fallback == fallback_before);
    CHECK(estimator.diagnostics().cost == cost_before);

    // The two channels side by side, which is the whole point of having two.
    // A step in which the EMBEDDED FILTER produced the estimate -- the fill-up
    // steps asserted above -- RETURNS SUCCESS and says so through the fallback
    // flag. A step in which NOTHING produced an estimate returns a failure. A
    // caller can act on the difference; before the conversion both arrived as
    // the same status enumerator on the same aggregate, separated only by an
    // undocumented conjunction with a flag that meant two different things.
    CHECK(fallback_before == false);
    CHECK_FALSE(refused.has_value());

    // The refusal does not latch, but the invalidated window must refill before
    // another optimization is allowed.
    estimator.predict(u);
    ctrlpp::Vector<double, NY> z_good;
    z_good << 0.7;
    REQUIRE(estimator.update(z_good).has_value());

    CHECK(estimator.diagnostics().status == ctrlpp::solve_status::optimal);
    CHECK(estimator.diagnostics().used_ekf_fallback);
    CHECK_FALSE(estimator.is_initialized());
}

TEST_CASE("MHE with inconsistent measurements", "[mhe][hardening][negative]")
{
    auto estimate_from = [](double sign) {
        ctrlpp::mhe_config<double, NX, NU, NY, N> cfg;
        cfg.Q = ctrlpp::Matrix<double, NX, NX>::Identity() * 0.01;
        cfg.R = ctrlpp::Matrix<double, NY, NY>::Identity() * 0.1;

        auto estimator = ctrlpp::test::constructed(MheType::create(linear_dynamics{}, position_measurement{}, cfg));
        ctrlpp::Vector<double, NU> u = ctrlpp::Vector<double, NU>::Zero();

        // Feed wildly inconsistent measurements (jumping between +100 and -100)
        for(std::size_t i = 0; i < 2 * N; ++i)
        {
            estimator.predict(u);
            ctrlpp::Vector<double, NY> z;
            z << sign * ((i % 2 == 0) ? 100.0 : -100.0);
            REQUIRE(estimator.update(z).has_value());
        }
        return estimator;
    };

    auto estimator = estimate_from(1.0);
    auto mirrored = estimate_from(-1.0);

    // Inconsistent data is not a failure and must not be reported as one: the
    // window problem is still convex and still solved, so the status is optimal
    // and the fallback is cleared. That disposition is the case's real subject
    // and no assertion on the estimate's magnitude can see it.
    CHECK(estimator.diagnostics().status == ctrlpp::solve_status::optimal);
    CHECK_FALSE(estimator.diagnostics().used_ekf_fallback);

    // The cost is quadratic with no constant offset and the dynamics are linear
    // with a zero input, so the whole estimation problem is ODD in the record:
    // negating every measurement negates the arrival prior, the residuals and
    // therefore the minimizer. Bitwise, because the negation of a double is
    // exact and every operation between the two runs is the same operation on
    // negated operands.
    CHECK(estimator.state() == -mirrored.state());

    // And the reported trajectory really is the estimate's own: its last node is
    // bitwise what `state()` returns, so a caller reading either gets the same
    // answer.
    const auto trajectory = estimator.trajectory();
    REQUIRE(trajectory.size() == N + 1);
    CHECK(trajectory[trajectory.size() - 1] == estimator.state());
}

// ── MHE hardening: precision ───────────────────────────────────────────────────

TEST_CASE("MHE linear system matches Kalman-like estimate", "[mhe][hardening][precision]")
{
    ctrlpp::mhe_config<double, NX, NU, NY, N> cfg;
    cfg.Q = ctrlpp::Matrix<double, NX, NX>::Identity() * 0.01;
    cfg.R = ctrlpp::Matrix<double, NY, NY>::Identity() * 0.1;
    cfg.P0 = ctrlpp::Matrix<double, NX, NX>::Identity() * 10.0;

    auto estimator = ctrlpp::test::constructed(MheType::create(linear_dynamics{}, position_measurement{}, cfg));

    // The equivalent linear filter the case name promises, built here rather
    // than assumed: the same transition and input matrices the dynamics model
    // differentiates to, the same measurement row, the same noise weights and
    // the same initial condition.
    ctrlpp::Matrix<double, NX, NX> A;
    A << 1.0, dt, 0.0, 1.0;
    ctrlpp::Matrix<double, NX, NU> B;
    B << 0.5 * dt * dt, dt;
    ctrlpp::Matrix<double, NY, NX> H;
    H << 1.0, 0.0;
    const ctrlpp::Matrix<double, NY, NU> D = ctrlpp::Matrix<double, NY, NU>::Zero();

    auto reference = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, NX, NU, NY>::create(
        ctrlpp::discrete_state_space<double, NX, NU, NY>{A, B, H, D},
        ctrlpp::kalman_config<double, NX, NU, NY>{.Q = cfg.Q, .R = cfg.R, .x0 = cfg.x0, .P0 = cfg.P0}));

    // True state: position 0, velocity 1 (constant velocity)
    ctrlpp::Vector<double, NX> x_true;
    x_true << 0.0, 1.0;
    ctrlpp::Vector<double, NU> u = ctrlpp::Vector<double, NU>::Zero();

    std::mt19937 gen(42);
    std::normal_distribution<double> noise(0.0, 0.1);

    // The two estimators are NOT bitwise equal and the reason is structural, so
    // it is pinned rather than hidden. The window problem uses the embedded
    // filter's covariance at the window head as its arrival cost, which is the
    // covariance conditioned on data strictly before the window rather than the
    // exact conditional one the filter carries forward; the quadratic program's
    // own convergence tolerance sets what remains. Both factors are stated here
    // and the realized disagreement is recorded against them.
    constexpr double solver_absolute_tolerance = 1e-3;
    constexpr double arrival_cost_amplification = 1.0 / 0.01; // one over the process weight
    const double budget = solver_absolute_tolerance / arrival_cost_amplification;

    double worst = 0.0;

    for(std::size_t i = 0; i < 3 * N; ++i)
    {
        x_true = linear_dynamics{}(x_true, u);
        ctrlpp::Vector<double, NY> z;
        z << x_true(0) + noise(gen);

        estimator.predict(u);
        REQUIRE(estimator.update(z).has_value());
        reference.predict(u);
        REQUIRE(reference.update(z).has_value());

        if(i >= N)
        {
            REQUIRE(estimator.diagnostics().status == ctrlpp::solve_status::optimal);
            REQUIRE_FALSE(estimator.diagnostics().used_ekf_fallback);
            worst = std::max(worst, (estimator.state() - reference.state()).cwiseAbs().maxCoeff());
        }
    }

    // Agreement with the constructed filter, at every window-solved step and in
    // both coordinates -- not with the truth, which the case name never claimed.
    CHECK(worst <= budget);
    // The agreement is not trivial: it is orders below the state itself, so the
    // comparison is doing work rather than comparing two numbers near zero.
    CHECK(worst < std::abs(reference.state()(0)));
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

    for(std::size_t i = 0; i < 5 * N; ++i)
    {
        x_true = linear_dynamics{}(x_true, u);
        ctrlpp::Vector<double, NY> z;
        z << x_true(0) + noise(gen);

        estimator.predict(u);
        REQUIRE(estimator.update(z).has_value());

        // The warm-up contract, asserted for the whole run rather than sampled:
        // the embedded filter supplies the estimate for exactly the first N
        // updates and the window solve supplies every one after.
        REQUIRE(estimator.diagnostics().used_ekf_fallback == (i < N));

        if(i < N)
            continue;

        // Convergence stated against what the estimator itself claims to know.
        // The old assertion counted how many of twenty-five steps happened to
        // improve, which is a statistic about one generator seed: a filter that
        // oscillates while slowly diverging can score more than half. This one
        // holds at EVERY step and in BOTH coordinates.
        REQUIRE(std::abs(estimator.state()(0) - x_true(0))
                <= reported_standard_errors * std::sqrt(estimator.covariance()(0, 0)));
        REQUIRE(std::abs(estimator.state()(1) - x_true(1))
                <= reported_standard_errors * std::sqrt(estimator.covariance()(1, 1)));
    }

    // The initial error of five was genuinely removed, and the reported
    // uncertainty came down with it: the position standard deviation started at
    // ten and ends far under one.
    CHECK(std::abs(estimator.state()(0) - x_true(0)) < 1.0);
    CHECK(estimator.covariance()(0, 0) < cfg.P0(0, 0));
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
    constexpr double measurement_step = 0.1;

    for(std::size_t i = 0; i < 2 * N; ++i)
    {
        estimator.predict(u);
        ctrlpp::Vector<double, NY> z;
        z << static_cast<double>(i) * measurement_step;
        REQUIRE(estimator.update(z).has_value());
    }

    const double last_measurement = static_cast<double>(2 * N - 1) * measurement_step;

    // The disposition, which is the case's real subject and which nothing read
    // before. A process weight ten decades below its neighbor makes the window
    // problem badly conditioned, and the backend says so: it reports an
    // inaccurate solve rather than an optimal one, and it does NOT fall back.
    // That is the honest report for this data and it is what a caller would act
    // on.
    CHECK(estimator.diagnostics().status == ctrlpp::solve_status::solved_inaccurate);
    CHECK_FALSE(estimator.diagnostics().used_ekf_fallback);

    // And the estimate is what the weighting DICTATES rather than merely finite.
    // The two process variances differ by ten decades, so their inverses -- the
    // penalties the window cost applies to a change in each coordinate -- differ
    // by the same ten decades. The velocity is therefore allowed to travel
    // across the window only that fraction of what the position travels, which
    // is a statement about the RATIO of the weights and needs no chosen number.
    const double weight_conditioning = cfg.Q(0, 0) / cfg.Q(1, 1);
    const double position_travel = std::abs(estimator.state()(0) - estimator.arrival_state()(0));
    CHECK(std::abs(estimator.state()(1) - estimator.arrival_state()(1))
          <= position_travel / weight_conditioning);

    // The position, whose weight is ordinary, tracks the measurement ramp. The
    // residual it is allowed is set by the measurement weight against the
    // position process weight: with a measurement variance of 0.1 and a process
    // variance of 1, the fit trades one against the other and the residual
    // cannot exceed the geometric mean of the two, which is the scale at which
    // the two terms of the cost balance.
    const double residual_scale = std::sqrt(cfg.R(0, 0) * cfg.Q(0, 0));
    CHECK(std::abs(estimator.state()(0) - last_measurement) <= residual_scale);
}
