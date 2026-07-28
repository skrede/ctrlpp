// What the oracles in this file decide.
//
// The dynamics and measurement these cases use are LINEAR, which is what makes
// the file's central claim checkable: for a linear model the extended filter is
// algebraically the linear filter, so an independently constructed
// ctrlpp::kalman_filter on the same data is an exact reference and not an
// approximate one. Every convergence and robustness case below is stated
// against that reference or against a closed form, never against finiteness and
// never against a fitted threshold.
//
// The equivalence is asserted twice, and the pair is the point:
//
//  * With ANALYTIC Jacobians the two filters execute the same arithmetic, so
//    they agree to a counted-operation budget at the scale of the state.
//  * With CENTRAL-DIFFERENCE Jacobians -- which is what this filter uses unless
//    the model supplies its own -- they do not. The difference quotient's step
//    is the cube root of epsilon, so cancellation in the numerator limits each
//    Jacobian entry to about epsilon to the two-thirds power, some five decades
//    coarser than the arithmetic around it. That cost is asserted rather than
//    hidden: the budget takes the ACTUAL Jacobian discrepancy, computed in the
//    test at the same operating point through the same helper the filter calls,
//    and carries it through the counted propagation into the state.
//
// What they deliberately do not decide. Nothing here asserts anything about
// nonlinear tracking; a genuinely nonlinear model has no closed-form reference,
// and inventing one would substitute a second implementation for an oracle.
// Nothing here validates the finite-difference step size either -- it is the
// filter's documented configuration parameter and the cases read its
// consequences rather than judging its value.

#include "hardening_helpers.h"

#include "ctrlpp/detail/numerical_diff.h"

#include "ctrlpp/estimation/ekf.h"
#include "ctrlpp/estimation/kalman.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;

namespace {

constexpr double ekf_eps = std::numeric_limits<double>::epsilon();

// Rounded operations along the longest chain through one filter step for a
// two-state, one-output problem, counted exactly as for the linear filter:
// five for the state propagation, seven for the covariance propagation, four
// for the innovation, seven for the innovation covariance, five for the gain
// solve, two for the state correction and thirteen for the Joseph update.
constexpr int ekf_step_ops = 43;

// Rounded operations carrying a Jacobian discrepancy into the state, counted
// along the chain it travels. A perturbation of F enters the propagated
// covariance twice, once through each factor of F P F'; the innovation
// covariance and the gain each divide by a quantity carrying that perturbation,
// doubling it again; and the correction multiplies the gain by the innovation.
// Eight bounds the accumulated amplification from above.
constexpr int ekf_jacobian_chain_ops = 8;

// The symmetric eigensolver's backward perturbation, bounded by the number of
// entries of a two-by-two times epsilon times the matrix norm (Weyl).
constexpr int ekf_eig_ops = 4;

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

/// @brief The same model, carrying the Jacobians it could have computed by hand.
///
/// The filter takes this path only when the model satisfies the differentiable
/// concepts, which requires BOTH state and input Jacobians on the dynamics; a
/// model supplying only one silently falls back to differencing, which is why
/// the static assertions below are here and not left implicit.
struct linear_dynamics_with_jacobians : linear_dynamics
{
    auto jacobian_x(const ctrlpp::Vector<double, 2>&, const ctrlpp::Vector<double, 1>&) const
        -> ctrlpp::Matrix<double, 2, 2>
    {
        ctrlpp::Matrix<double, 2, 2> F;
        F << 1.0, 0.1, 0.0, 1.0;
        return F;
    }

    auto jacobian_u(const ctrlpp::Vector<double, 2>&, const ctrlpp::Vector<double, 1>&) const
        -> ctrlpp::Matrix<double, 2, 1>
    {
        ctrlpp::Matrix<double, 2, 1> G;
        G << 0.0, 0.1;
        return G;
    }
};

struct linear_measurement_with_jacobian : linear_measurement
{
    auto jacobian(const ctrlpp::Vector<double, 2>&) const -> ctrlpp::Matrix<double, 1, 2>
    {
        ctrlpp::Matrix<double, 1, 2> H;
        H << 1.0, 0.0;
        return H;
    }
};

static_assert(ctrlpp::differentiable_dynamics<linear_dynamics_with_jacobians, double, 2, 1>);
static_assert(ctrlpp::differentiable_measurement<linear_measurement_with_jacobian, double, 2, 1>);
static_assert(!ctrlpp::differentiable_dynamics<linear_dynamics, double, 2, 1>);

using ekf_t = ctrlpp::ekf<double, 2, 1, 1, linear_dynamics, linear_measurement>;
using ekf_analytic_t =
    ctrlpp::ekf<double, 2, 1, 1, linear_dynamics_with_jacobians, linear_measurement_with_jacobian>;

/// @brief The linear system the model above IS, written out so an independent
/// linear filter can be driven with the same data.
auto equivalent_system() -> ctrlpp::discrete_state_space<double, 2, 1, 1>
{
    ctrlpp::discrete_state_space<double, 2, 1, 1> sys;
    sys.A << 1.0, 0.1, 0.0, 1.0;
    sys.B << 0.0, 0.1;
    sys.C << 1.0, 0.0;
    sys.D << 0.0;
    return sys;
}

/// @brief Spectral norm of the closed-loop error map the reference filter's own
/// reported covariance implies, used to carry a divergence bound forward.
///
/// It exceeds one during the transient, so a bound formed as a plain multiple of
/// the step count would not be a bound at all.
auto closed_loop_norm(const ctrlpp::discrete_state_space<double, 2, 1, 1>& sys,
                      const Eigen::Matrix2d& P_pred, double R) -> double
{
    const double S = (sys.C * P_pred * sys.C.transpose())(0, 0) + R;
    const Eigen::Vector2d K = (P_pred * sys.C.transpose()).eval() / S;
    return ((Eigen::Matrix2d::Identity() - K * sys.C) * sys.A).norm();
}

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

TEST_CASE("EKF for a linear system IS the Kalman filter", "[ekf][hardening][precision]")
{
    // The name this case used to carry promised agreement with the Kalman filter
    // to one percent, and then constructed no Kalman filter: it asserted that the
    // position estimate was within 1.0 of the truth, which is one measurement
    // standard deviation written as a bare literal and true of a filter that
    // ignored its measurements entirely.
    const auto sys = equivalent_system();

    ctrlpp::ekf_config<double, 2, 1, 1> cfg{};
    cfg.Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    cfg.R = Eigen::Matrix<double, 1, 1>::Identity();
    cfg.x0 = Eigen::Vector2d::Zero();
    cfg.P0 = Eigen::Matrix<double, 2, 2>::Identity() * 10.0;

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;

    SECTION("with analytic Jacobians, to the shared arithmetic's own rounding")
    {
        auto filter = ctrlpp::test::constructed(ekf_analytic_t::create(
            linear_dynamics_with_jacobians{}, linear_measurement_with_jacobian{}, cfg));
        auto reference = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(
            sys, {.Q = cfg.Q, .R = cfg.R, .x0 = cfg.x0, .P0 = cfg.P0}));

        double true_pos = 0.0;
        double divergence_budget = 0.0;

        for(int k = 0; k < 50; ++k)
        {
            CAPTURE(k);
            true_pos += 0.1;
            filter.predict(u);
            reference.predict(u);

            const double amplification = closed_loop_norm(sys, reference.covariance(), cfg.R(0, 0));

            Eigen::Matrix<double, 1, 1> z;
            z << true_pos;
            REQUIRE(filter.update(z).has_value());
            REQUIRE(reference.update(z).has_value());

            const double scale = std::max({filter.state().cwiseAbs().maxCoeff(),
                                           reference.state().cwiseAbs().maxCoeff(), 1.0});
            divergence_budget =
                amplification * divergence_budget + ekf_step_ops * ekf_eps * scale;

            CAPTURE(filter.state()[0], reference.state()[0]);
            REQUIRE_THAT(filter.state()[0], WithinAbs(reference.state()[0], divergence_budget));
            REQUIRE_THAT(filter.state()[1], WithinAbs(reference.state()[1], divergence_budget));
        }
    }

    SECTION("with central-difference Jacobians, to the differencing error they carry")
    {
        auto filter = ctrlpp::test::constructed(
            ekf_t::create(linear_dynamics{}, linear_measurement{}, cfg));
        auto reference = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(
            sys, {.Q = cfg.Q, .R = cfg.R, .x0 = cfg.x0, .P0 = cfg.P0}));

        double true_pos = 0.0;
        double divergence_budget = 0.0;
        double worst_jacobian_error = 0.0;

        for(int k = 0; k < 50; ++k)
        {
            CAPTURE(k);
            true_pos += 0.1;

            // The discrepancy the filter is about to introduce, measured at the
            // state it is about to linearize around, through the same helper it
            // calls and with the same step it uses by default. This is what the
            // budget is built from, so the budget describes the filter's actual
            // approximation rather than a worst case guessed from the outside.
            const auto F = ctrlpp::detail::numerical_jacobian_x<double, 2, 1>(
                linear_dynamics{}, filter.state(), u, cfg.numerical_eps);
            const auto H = ctrlpp::detail::numerical_jacobian_h<double, 2, 1>(
                linear_measurement{}, filter.state(), cfg.numerical_eps);
            const double jacobian_error =
                (F - sys.A).cwiseAbs().maxCoeff() + (H - sys.C).cwiseAbs().maxCoeff();
            worst_jacobian_error = std::max(worst_jacobian_error, jacobian_error);

            filter.predict(u);
            reference.predict(u);

            const double amplification = closed_loop_norm(sys, reference.covariance(), cfg.R(0, 0));

            Eigen::Matrix<double, 1, 1> z;
            z << true_pos;
            REQUIRE(filter.update(z).has_value());
            REQUIRE(reference.update(z).has_value());

            const double scale = std::max({filter.state().cwiseAbs().maxCoeff(),
                                           reference.state().cwiseAbs().maxCoeff(), 1.0});
            // The floor at epsilon covers the steps where the difference quotient
            // happens to be exact -- at the first step it is, because the state is
            // exactly zero -- and the arithmetic rounding still applies there.
            divergence_budget = amplification * divergence_budget
                                + ekf_jacobian_chain_ops * std::max(jacobian_error, ekf_eps) * scale;

            CAPTURE(jacobian_error, filter.state()[0], reference.state()[0]);
            REQUIRE_THAT(filter.state()[0], WithinAbs(reference.state()[0], divergence_budget));
            REQUIRE_THAT(filter.state()[1], WithinAbs(reference.state()[1], divergence_budget));
        }

        // The measurement model is the first state coordinate, so its difference
        // quotient is (x+h) - (x-h) over 2h: EXACT, with no cancellation at all.
        // The whole of the discrepancy therefore comes from the dynamics, and it
        // sits at the cube-root-of-epsilon scale that the step size fixes.
        CAPTURE(worst_jacobian_error);
        REQUIRE(worst_jacobian_error > ekf_eps);
        REQUIRE(worst_jacobian_error < ekf_jacobian_chain_ops * std::pow(ekf_eps, 2.0 / 3.0));
    }
}

TEST_CASE("EKF covariance stays PD over 1000 steps", "[ekf][hardening][stability]")
{
    auto filter = make_ekf();

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
        // propagated and the corrected covariance, so the two off-diagonal
        // entries are the same double.
        REQUIRE(filter.covariance()(0, 1) == filter.covariance()(1, 0));

        // The floor is the eigensolver's own backward error at the scale of the
        // matrix handed to it. The threshold it replaces, -1e-10, is about
        // 450000 epsilons of slack BELOW zero, so a genuinely indefinite
        // covariance passed for a thousand consecutive steps.
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> eigsolver(filter.covariance());
        REQUIRE(eigsolver.eigenvalues().minCoeff()
                >= -ekf_eig_ops * ekf_eps * filter.covariance().norm());
    }
}

TEST_CASE("EKF estimation error follows its own closed-loop recursion",
          "[ekf][hardening][convergence]")
{
    const auto sys = equivalent_system();
    auto filter = make_ekf();

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;

    // The measurements are the observed coordinate of a state advanced by the
    // model itself, with no noise, so the estimation error is homogeneous:
    //     e(k+1) = (I - K(k) C) A e(k).
    // The test propagates it alongside using the gain implied by the filter's own
    // reported predicted covariance, and the oracle is that propagated error --
    // not the fitted 0.5 it replaces, which sat eight decades above the realized
    // value and would have accepted a filter that never converged at all.
    Eigen::Vector2d x_true;
    x_true << 0.0, 1.0;
    Eigen::Vector2d error = Eigen::Vector2d::Zero() - x_true;

    // The Jacobian discrepancy enters this recursion the same way it enters the
    // equivalence case, so the budget is built the same way.
    double divergence_budget = 0.0;

    for(int k = 0; k < 200; ++k)
    {
        CAPTURE(k);
        // make_ekf leaves the differencing step at its default, which is the
        // cube root of epsilon, so the discrepancy is measured with that step.
        const auto F = ctrlpp::detail::numerical_jacobian_x<double, 2, 1>(
            linear_dynamics{}, filter.state(), u, std::cbrt(ekf_eps));
        const double jacobian_error = (F - sys.A).cwiseAbs().maxCoeff();

        x_true = (sys.A * x_true).eval();
        filter.predict(u);

        const Eigen::Matrix2d P_pred = filter.covariance();
        const double S = (sys.C * P_pred * sys.C.transpose())(0, 0) + 1.0;
        const Eigen::Vector2d K = (P_pred * sys.C.transpose()).eval() / S;
        const Eigen::Matrix2d closed_loop = ((Eigen::Matrix2d::Identity() - K * sys.C) * sys.A).eval();
        error = (closed_loop * error).eval();

        Eigen::Matrix<double, 1, 1> z;
        z << x_true(0);
        REQUIRE(filter.update(z).has_value());

        const double scale = std::max({x_true.cwiseAbs().maxCoeff(),
                                       filter.state().cwiseAbs().maxCoeff(), 1.0});
        divergence_budget = closed_loop.norm() * divergence_budget
                            + ekf_jacobian_chain_ops * std::max(jacobian_error, ekf_eps) * scale;

        const Eigen::Vector2d realized = filter.state() - x_true;
        CAPTURE(realized(0), realized(1), error(0), error(1));
        REQUIRE_THAT(realized(0), WithinAbs(error(0), divergence_budget));
        REQUIRE_THAT(realized(1), WithinAbs(error(1), divergence_budget));
    }
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

    // The reference is the linear filter on the identical ill-conditioned data.
    // Conditioning is exactly the thing finiteness could not see: a filter whose
    // near-singular process noise had destroyed its gain would still report
    // finite numbers forever, and would immediately part company with the
    // reference.
    const auto sys = equivalent_system();
    auto reference = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(
        sys, {.Q = cfg.Q, .R = cfg.R, .x0 = cfg.x0, .P0 = cfg.P0}));

    Eigen::Matrix<double, 1, 1> u;
    u << 0.0;
    double divergence_budget = 0.0;

    for(int k = 0; k < 100; ++k)
    {
        CAPTURE(k);
        const auto F = ctrlpp::detail::numerical_jacobian_x<double, 2, 1>(
            linear_dynamics{}, filter.state(), u, cfg.numerical_eps);
        const double jacobian_error = (F - sys.A).cwiseAbs().maxCoeff();

        filter.predict(u);
        reference.predict(u);

        const double amplification = closed_loop_norm(sys, reference.covariance(), cfg.R(0, 0));

        Eigen::Matrix<double, 1, 1> z;
        z << 1.0;
        REQUIRE(filter.update(z).has_value());
        REQUIRE(reference.update(z).has_value());

        const double scale = std::max({filter.state().cwiseAbs().maxCoeff(),
                                       reference.state().cwiseAbs().maxCoeff(), 1.0});
        divergence_budget = amplification * divergence_budget
                            + ekf_jacobian_chain_ops * std::max(jacobian_error, ekf_eps) * scale;

        CAPTURE(filter.state()[0], reference.state()[0]);
        REQUIRE_THAT(filter.state()[0], WithinAbs(reference.state()[0], divergence_budget));
        REQUIRE_THAT(filter.state()[1], WithinAbs(reference.state()[1], divergence_budget));

        REQUIRE(filter.covariance()(0, 1) == filter.covariance()(1, 0));
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> eigsolver(filter.covariance());
        REQUIRE(eigsolver.eigenvalues().minCoeff()
                >= -ekf_eig_ops * ekf_eps * filter.covariance().norm());
    }

    // The measurement is constant, so the observed coordinate must have settled
    // on it to within the uncertainty the filter itself reports. A multiplier of
    // one is the consistency contract, not a fitted margin: an error larger than
    // the reported standard deviation means the filter is claiming a confidence
    // its own estimate does not earn, which is the defect this case exists to
    // catch.
    CAPTURE(filter.state()[0], filter.covariance()(0, 0));
    REQUIRE(std::abs(filter.state()[0] - 1.0) <= std::sqrt(filter.covariance()(0, 0)));
}
