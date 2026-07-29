// What the oracles in this file decide.
//
// The observer carries a state and no covariance, so every claim here is a
// claim about that state and about the typed rejection that protects it:
//
//  * A rejected measurement leaves the carried state BITWISE unchanged, and a
//    later valid step produces exactly what it would have produced had the bad
//    sample never arrived. Both are exact comparisons, because a rejected step
//    performs no arithmetic at all and a tolerance would admit one that
//    partially ran.
//  * With a zero gain the observer IS the open-loop system, so its state after
//    ten steps is the tenth power of the system matrix applied to the initial
//    state -- one coordinate bitwise unchanged, the other within a counted
//    budget.
//  * Convergence is asserted against the observer's OWN closed-loop matrix
//    (I - L C) A, raised to the number of steps taken and applied to the initial
//    error. Because the true state of these cases is exactly zero and stays
//    exactly zero, the observer's state IS the error, and the recursion for it
//    is exact. The oracle is therefore the realized error itself, not a decay
//    rate and not a fitted final threshold.
//
// What they deliberately do not decide. Nothing here asserts a pole placement:
// the gain is given, not designed, and the cases read the closed-loop matrix
// that gain produces rather than checking that it is a good one. Nothing here
// asserts anything about the observer under a non-zero input either; every
// convergence case drives it with a zero input so that the error recursion is
// homogeneous and the closed form applies without an inhomogeneous term.

#include "hardening_helpers.h"
#include "ctrlpp/estimation/luenberger.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

constexpr double obs_eps = std::numeric_limits<double>::epsilon();

// Rounded operations in one observer step, counted along the longest chain
// reaching one state component. The prediction contracts a two-column row of A
// against the state (two multiplies and one addition) and adds the input term
// (one multiply and one addition): five. The correction contracts C against the
// predicted state (two multiplies and one addition), subtracts it from the
// measurement (one), scales by the gain component (one) and adds it back (one):
// five more.
constexpr int observer_step_ops = 10;

auto make_system()
{
    ctrlpp::discrete_state_space<double, 2, 1, 1> sys;
    sys.A << 0.9, 0.1, 0.0, 0.8;
    sys.B << 0.0, 1.0;
    sys.C << 1.0, 0.0;
    sys.D << 0.0;
    return sys;
}

/// @brief The observer's own closed-loop error map, (I - L C) A.
///
/// Derived rather than assumed. With the loop ordering these cases use --
/// predict, then correct with a measurement of the CURRENT true state, then
/// advance the truth -- the error obeys
///   e(k+1) = (I - L C) A e(k) - (A - L C) x_true(k),
/// and every case below drives the truth from exactly zero with a zero input, so
/// the inhomogeneous term vanishes identically and the error recursion is
/// exactly e(k+1) = (I - L C) A e(k).
auto closed_loop_map(const ctrlpp::discrete_state_space<double, 2, 1, 1>& sys,
                     const Eigen::Matrix<double, 2, 1>& L) -> Eigen::Matrix2d
{
    return ((Eigen::Matrix2d::Identity() - L * sys.C) * sys.A).eval();
}

auto matrix_power(Eigen::Matrix2d M, int n) -> Eigen::Matrix2d
{
    Eigen::Matrix2d result = Eigen::Matrix2d::Identity();
    for(int k = 0; k < n; ++k)
        result = (result * M).eval();
    return result;
}

}

TEST_CASE("Luenberger NaN measurement is rejected without touching the state",
          "[luenberger][hardening][negative]")
{
    auto sys = make_system();
    Eigen::Matrix<double, 2, 1> L;
    // Both gain entries are nonzero, so a measurement that reached the fused
    // correction x + L*(z - C*x) would poison BOTH state components. The
    // observer carries no covariance, so the preserved-invariant argument has a
    // state half only.
    L << 0.5, 0.3;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();

    ctrlpp::luenberger_observer<double, 2, 1, 1> obs(sys, L, x0);
    // Stepped only with the valid measurement, never with the poisoned one, so
    // it says what the observer would have carried had the bad sample never
    // arrived.
    ctrlpp::luenberger_observer<double, 2, 1, 1> reference(sys, L, x0);

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;
    obs.predict(u);
    reference.predict(u);

    // Snapshot immediately before the poisoned step.
    const Eigen::Vector2d x_before = obs.state();

    Eigen::Matrix<double, 1, 1> z_bad;
    z_bad << std::numeric_limits<double>::quiet_NaN();

    const auto rejected = obs.update(z_bad);

    REQUIRE_FALSE(rejected.has_value());
    REQUIRE(rejected.error() == ctrlpp::luenberger_update_error::non_finite_measurement);

    // Exact comparison, not a tolerance: a rejected step performs no arithmetic
    // on the carried state at all, so bitwise equality is the contract and a
    // tolerance would admit a step that partially ran.
    CHECK(obs.state() == x_before);
    // A rejection describes the sample, not the observer: nothing was mutated,
    // so the observer is not degraded and must not report that it is.
    CHECK(obs.health() == ctrlpp::luenberger_health::ok);

    // The poison did not latch: the next valid step produces exactly what it
    // would have produced had the poisoned step never been attempted.
    Eigen::Matrix<double, 1, 1> z_good;
    z_good << 1.0;
    REQUIRE(obs.update(z_good).has_value());
    REQUIRE(reference.update(z_good).has_value());

    CHECK(obs.state() == reference.state());
}

TEST_CASE("Luenberger rejects correction overflow without committing it",
          "[luenberger][hardening][negative]")
{
    auto sys = make_system();
    Eigen::Matrix<double, 2, 1> gain;
    gain << std::numeric_limits<double>::max(), 0.0;
    ctrlpp::luenberger_observer<double, 2, 1, 1> observer{
        sys, gain, Eigen::Vector2d::Zero()};
    auto const state_before = observer.state();

    Eigen::Matrix<double, 1, 1> measurement;
    measurement << 2.0;
    auto result = observer.update(measurement);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error()
          == ctrlpp::luenberger_update_error::non_finite_result);
    CHECK(observer.state() == state_before);
    CHECK(observer.health() == ctrlpp::luenberger_health::ok);
}

TEST_CASE("Luenberger zero observer gains (open-loop)",
          "[luenberger][hardening][negative]")
{
    auto sys = make_system();
    Eigen::Matrix<double, 2, 1> L = Eigen::Matrix<double, 2, 1>::Zero();
    Eigen::Vector2d x0;
    x0 << 1.0, 0.0;

    ctrlpp::luenberger_observer<double, 2, 1, 1> obs(sys, L, x0);

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;

    // With L=0, update has no effect -- purely open loop
    constexpr int steps = 10;
    for(int k = 0; k < steps; ++k)
    {
        obs.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << 5.0;
        REQUIRE(obs.update(z).has_value());
    }

    // A zero gain means no correction, so the observer is the open-loop system
    // and its state is A^10 applied to the initial state. That is what actually
    // proves the gain is doing nothing; finiteness held for any correction at
    // all, including a fully engaged one, so it could not see the property the
    // case is named for. The measurement of 5.0 is deliberately far from the
    // trajectory: an observer that applied even a fraction of it would move
    // toward five and fail both assertions below.

    // Bitwise: the second row of A is [0, 0.8] and the second component starts
    // at exactly zero, so every step multiplies zero by 0.8 and adds an exact
    // zero product.
    REQUIRE(obs.state()[1] == 0.0);

    // The first component is 0.9^10 reached by ten multiplications -- the
    // coupling term 0.1 * x1 is an exact zero at every step -- compared against
    // one call to the standard power function, whose own error is bounded by one
    // ulp. Eleven rounded operations in all.
    constexpr int open_loop_ops = steps + 1;
    REQUIRE_THAT(obs.state()[0], WithinRel(std::pow(0.9, steps), open_loop_ops * obs_eps));
}

TEST_CASE("Luenberger observer error decays for stable poles",
          "[luenberger][hardening][stability]")
{
    auto sys = make_system();
    // Choose L to place observer poles well inside unit circle
    Eigen::Matrix<double, 2, 1> L;
    L << 0.5, 0.3;
    Eigen::Vector2d x0;
    x0 << 10.0, 5.0;

    ctrlpp::luenberger_observer<double, 2, 1, 1> obs(sys, L, x0);

    // True state
    Eigen::Vector2d x_true;
    x_true << 0.0, 0.0;
    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;

    constexpr int steps = 200;
    for(int k = 0; k < steps; ++k)
    {
        obs.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << sys.C(0, 0) * x_true[0] + sys.C(0, 1) * x_true[1];
        REQUIRE(obs.update(z).has_value());

        x_true = (sys.A * x_true + sys.B * u).eval();
    }

    // The truth started at exactly zero and is driven by a zero input through a
    // matrix product, so it is still exactly zero and the error recursion is
    // homogeneous. Assert that first, because the closed form below depends on it.
    REQUIRE(x_true == Eigen::Vector2d::Zero());

    // The error after 200 steps is the closed-loop map raised to 200 and applied
    // to the initial error -- not merely "smaller than it started" (which the
    // old comparison against the INITIAL error, misnamed prev_error, asserted)
    // and not "below 0.1" (which was a fitted constant three decades above the
    // realized value). The relative budget is 200 steps of the observer's own
    // arithmetic plus the same number of multiplications forming the reference
    // power: a perturbation injected at step k is itself contracted by the
    // remaining 200-k applications of the map, so the relative error accumulates
    // additively rather than being swamped by the initial scale.
    const Eigen::Matrix2d M = closed_loop_map(sys, L);
    const Eigen::Vector2d expected_error = matrix_power(M, steps) * x0;
    const double decay_budget = 2 * steps * observer_step_ops * obs_eps;

    const Eigen::Vector2d error = obs.state() - x_true;
    CAPTURE(error(0), error(1), expected_error(0), expected_error(1));
    REQUIRE_THAT(error(0), WithinRel(expected_error(0), decay_budget));
    REQUIRE_THAT(error(1), WithinRel(expected_error(1), decay_budget));
}

TEST_CASE("Luenberger observer converges to true state within 200 steps",
          "[luenberger][hardening][convergence]")
{
    auto sys = make_system();
    Eigen::Matrix<double, 2, 1> L;
    L << 0.5, 0.3;
    Eigen::Vector2d x0;
    x0 << 5.0, 2.0;

    ctrlpp::luenberger_observer<double, 2, 1, 1> obs(sys, L, x0);

    Eigen::Vector2d x_true;
    x_true << 0.0, 0.0;
    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;

    constexpr int steps = 200;
    for(int k = 0; k < steps; ++k)
    {
        obs.predict(u);
        double meas = sys.C(0, 0) * x_true[0] + sys.C(0, 1) * x_true[1];
        Eigen::Matrix<double, 1, 1> z;
        z << meas;
        REQUIRE(obs.update(z).has_value());

        x_true = (sys.A * x_true + sys.B * u).eval();
    }

    REQUIRE(x_true == Eigen::Vector2d::Zero());

    // Same closed form, a different initial error. The threshold this replaces
    // was 0.01 against a realized error of order 1e-29 -- twenty-seven decades of
    // slack, which is a bound only in name. Asserting the closed form instead
    // fails an observer that converges too SLOWLY and equally one that converges
    // faster than its own closed-loop matrix permits, which no one-sided
    // threshold can do.
    const Eigen::Matrix2d M = closed_loop_map(sys, L);
    const Eigen::Vector2d expected_error = matrix_power(M, steps) * x0;
    const double decay_budget = 2 * steps * observer_step_ops * obs_eps;

    CAPTURE(obs.state()[0], obs.state()[1], expected_error(0), expected_error(1));
    REQUIRE_THAT(obs.state()[0] - x_true[0], WithinRel(expected_error(0), decay_budget));
    REQUIRE_THAT(obs.state()[1] - x_true[1], WithinRel(expected_error(1), decay_budget));
}
