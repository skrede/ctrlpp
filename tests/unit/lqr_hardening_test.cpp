// What the oracles in this file decide.
//
// The gain helper forwards the Riccati solver's verdict, so every claim here is
// a claim about that verdict or about the gain the solution implies:
//
//  * A weighting the solver cannot accept is refused with the enumerator that
//    names the cause, never with a gain the caller could mistake for usable.
//  * A well-posed pair is SOLVED, and the solution is held to the equation it
//    is supposed to satisfy -- the Riccati residual against a counted-operation
//    budget scaled by the largest of the four terms that cancel to produce it --
//    not merely to being finite.
//  * The returned gain is the one that solution implies, formed independently in
//    the test as (R + B'PB)^-1 B'PA, and the closed loop it produces is
//    asymptotically stable. Strictly inside the unit circle is the exact
//    contract boundary, not a fitted constant.
//  * Where the problem has a closed form (the scalar integrator's golden-ratio
//    fixed point) the closed form is asserted, with the budget carried from the
//    residual through the residual's own Frechet derivative at the solution.
//  * The sequence helpers refuse every length, step or entry they cannot index
//    or compute on, with the enumerator that names it. This file is built with
//    the standard library's bounds assertions, so an out-of-range read that a
//    refusal misses aborts the test instead of returning garbage.
//
// What they deliberately do not decide. Nothing here asserts optimality against
// a competing controller, and nothing asserts a conditioning model: the
// ill-conditioned weighting case in the sibling Riccati file records why a
// counted-operation budget is the wrong oracle in that regime rather than
// widening one until it passes.

#include "hardening_helpers.h"
#include "ctrlpp/control/lqr.h"
#include "ctrlpp/control/dare.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>
#include <span>
#include <vector>
#include <cstddef>
#include <utility>

namespace
{

using state_matrix = Eigen::Matrix<double, 2, 2>;
using input_matrix = Eigen::Matrix<double, 2, 1>;
using weight_matrix = Eigen::Matrix<double, 1, 1>;
using gain_matrix = Eigen::Matrix<double, 1, 2>;
using tv_result = ctrlpp::expected<std::vector<gain_matrix>, ctrlpp::lqr_error>;

struct tv_sequences
{
    std::vector<state_matrix> As;
    std::vector<input_matrix> Bs;
    std::vector<state_matrix> Qs;
    std::vector<weight_matrix> Rs;
    state_matrix Qf;
};

// The double integrator with unit weights, each per-step sequence sized on its own.
auto double_integrator_sequences(std::size_t na, std::size_t nb, std::size_t nq, std::size_t nr) -> tv_sequences
{
    state_matrix A;
    A << 1.0, 1.0, 0.0, 1.0;
    input_matrix B;
    B << 0.5, 1.0;
    const weight_matrix R = weight_matrix::Identity();
    const state_matrix Q = state_matrix::Identity();
    return {std::vector<state_matrix>(na, A), std::vector<input_matrix>(nb, B), std::vector<state_matrix>(nq, Q),
            std::vector<weight_matrix>(nr, R), Q};
}

auto tv_gains_of(const tv_sequences& s) -> tv_result
{
    return ctrlpp::lqr_tv_gains<double, 2, 1>(s.As, s.Bs, s.Qs, s.Rs, s.Qf);
}

auto cost_of(const std::vector<input_matrix>& xs, const std::vector<weight_matrix>& us, const state_matrix& Q,
             const weight_matrix& R) -> ctrlpp::expected<double, ctrlpp::lqr_error>
{
    return ctrlpp::lqr_cost<double, 2, 1>(std::span<const input_matrix>{xs}, std::span<const weight_matrix>{us}, Q, R);
}

// Dyadic entries small enough that every product and partial sum is exact in
// double, so the trajectory cost has a single correct value and no tolerance.
auto dyadic_states(std::size_t n) -> std::vector<input_matrix>
{
    std::vector<input_matrix> xs(n);
    for(std::size_t k = 0; k < n; ++k)
        xs[k] << 0.5 * static_cast<double>(k + 1), -0.25 * static_cast<double>(k);
    return xs;
}

auto dyadic_inputs(std::size_t n) -> std::vector<weight_matrix>
{
    std::vector<weight_matrix> us(n);
    for(std::size_t k = 0; k < n; ++k)
        us[k] << 0.125 * static_cast<double>(k + 2);
    return us;
}

auto trajectory_weights() -> std::pair<state_matrix, weight_matrix>
{
    return {state_matrix::Identity(), 2.0 * weight_matrix::Identity()};
}

auto ramp_gains(std::size_t n) -> std::vector<gain_matrix>
{
    std::vector<gain_matrix> gains(n);
    for(std::size_t k = 0; k < n; ++k)
        gains[k] << 1.0 + static_cast<double>(k), -0.5 * static_cast<double>(k);
    return gains;
}

}

TEST_CASE("LQR refuses a non-finite state weighting", "[lqr][hardening][negative]")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, 1.0, 0.0, 1.0;
    B << 0.5, 1.0;
    Q = ctrlpp::test::nan_matrix<double, 2, 2>();
    R << 1.0;

    auto result = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    // The case name already claimed the answer; this holds the library to it.
    // A gain formed from a NaN weighting would be a controller reporting
    // success while carrying no usable number, so "nullopt or NaN" was a
    // disjunction over the only two outcomes observable here and could not
    // fail whichever the library did.
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::dare_error::non_finite_input);
}

TEST_CASE("LQR refuses a non-finite input weighting", "[lqr][hardening][negative]")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, 1.0, 0.0, 1.0;
    B << 0.5, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R << std::numeric_limits<double>::quiet_NaN();

    auto result = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::dare_error::non_finite_input);
}

TEST_CASE("LQR stabilizes an ill-conditioned but well-posed pair", "[lqr][hardening][robustness]")
{
    // A = diag(1, 1e-10), B = [1; 0]. Mode 0 sits at eigenvalue 1 and IS input
    // coupled, so it is controllable; mode 1 sits at 1e-10, which is inside the
    // unit circle and so needs no control at all. The pair is therefore
    // stabilizable (though not controllable -- [B, AB] has rank 1), and Q = I
    // makes it detectable, so a unique stabilizing solution exists and a
    // refusal here would be a solvable problem reported unsolvable.
    auto A = ctrlpp::test::ill_conditioned_2x2<double>(1e10);
    Eigen::Matrix<double, 2, 1> B;
    B << 1.0, 0.0;
    auto Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto const solved = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(solved.has_value());

    auto const& P = solved->P;
    constexpr double eps = std::numeric_limits<double>::epsilon();

    // The solution solves the equation, not merely "is finite".
    auto const res = ctrlpp::test::riccati_residual<double, 2, 1>(A, B, Q, R, P);
    CAPTURE(res.norm, res.scale);
    REQUIRE(res.norm <= ctrlpp::test::riccati_residual_ops<2, 1> * eps * res.scale);

    // ... and it is positive definite.
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> pes(P);
    for(int i = 0; i < 2; ++i)
        CHECK(pes.eigenvalues()(i) > 0.0);

    auto const result = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    // The gain is the one this solution implies, formed the same way.
    auto const K_expected = ctrlpp::test::riccati_gain<double, 2, 1>(A, B, R, P);
    CHECK((*result - K_expected).norm()
          <= ctrlpp::test::riccati_residual_ops<2, 1> * eps * K_expected.norm());

    // The mathematically required property of an LQR gain: the closed loop is
    // asymptotically stable. Strictly inside the unit circle is the exact
    // contract boundary, not a fitted constant.
    Eigen::Matrix<double, 2, 2> Acl = (A - B * (*result)).eval();
    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> ces(Acl, false);
    for(int i = 0; i < 2; ++i)
        REQUIRE(std::abs(ces.eigenvalues()(i)) < 1.0);
}

TEST_CASE("LQR double integrator matches analytical", "[lqr][hardening][precision]")
{
    double dt = 0.1;
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, dt, 0.0, 1.0;
    B << 0.5 * dt * dt, dt;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R << 1.0;

    auto result = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    auto const& K = *result;
    constexpr double eps = std::numeric_limits<double>::epsilon();

    // The case name promises the analytical gain and nothing analytical was
    // asserted: finiteness held for every gain the helper could conceivably
    // return, including one that destabilizes the loop. The gain is fixed by the
    // Riccati solution for this data, so the solution is checked against the
    // equation first and the gain against the solution second.
    auto const solved = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(solved.has_value());

    auto const res = ctrlpp::test::riccati_residual<double, 2, 1>(A, B, Q, R, solved->P);
    CAPTURE(res.norm, res.scale);
    REQUIRE(res.norm <= ctrlpp::test::riccati_residual_ops<2, 1> * eps * res.scale);

    auto const K_expected = ctrlpp::test::riccati_gain<double, 2, 1>(A, B, R, solved->P);
    CAPTURE(K(0, 0), K(0, 1), K_expected(0, 0), K_expected(0, 1));
    CHECK((K - K_expected).norm()
          <= ctrlpp::test::riccati_residual_ops<2, 1> * eps * K_expected.norm());

    // Verify closed-loop eigenvalues are inside unit circle
    auto Acl = (A - B * K).eval();
    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> solver(Acl, false);
    for(int i = 0; i < 2; ++i)
        CHECK(std::abs(solver.eigenvalues()(i)) < 1.0);
}

TEST_CASE("LQR closed-loop eigenvalues inside unit circle", "[lqr][hardening][stability]")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, 1.0, 0.0, 1.0;
    B << 0.5, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R << 1.0;

    auto result = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    auto Acl = (A - B * *result).eval();
    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> solver(Acl, false);
    for(int i = 0; i < 2; ++i)
        REQUIRE(std::abs(solver.eigenvalues()(i)) < 1.0);

    // Stability alone does not say the gain is the OPTIMAL one: a detuned gain
    // that still places both eigenvalues inside the circle passes the loop
    // above. Optimality is the Riccati equation, so it is asserted here too.
    constexpr double eps = std::numeric_limits<double>::epsilon();
    auto const solved = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(solved.has_value());

    auto const res = ctrlpp::test::riccati_residual<double, 2, 1>(A, B, Q, R, solved->P);
    CAPTURE(res.norm, res.scale);
    REQUIRE(res.norm <= ctrlpp::test::riccati_residual_ops<2, 1> * eps * res.scale);

    auto const K_expected = ctrlpp::test::riccati_gain<double, 2, 1>(A, B, R, solved->P);
    CHECK((*result - K_expected).norm()
          <= ctrlpp::test::riccati_residual_ops<2, 1> * eps * K_expected.norm());
}

TEST_CASE("LQR refuses an unstabilizable pair with the enumerator that names it",
          "[lqr][hardening][negative]")
{
    // A = diag(1e-10, 1), B = [1; 0]. Mode 1 sits at eigenvalue exactly 1 with
    // B(1) = 0: uncontrollable AND not asymptotically stable. The pair is
    // therefore NOT stabilizable and no stabilizing solution exists, so a gain
    // returned here would be a controller claiming to stabilize a plant it
    // cannot.
    //
    // The enumerator is reachable in this case because the uncontrollable mode
    // is ON the unit circle: it contributes 1 and its reciprocal 1 to the
    // symplectic spectrum, neither strictly inside, so only one of the four
    // eigenvalues lands in the stable region and the placement count falls
    // short of the state dimension.
    auto sys = ctrlpp::test::near_singular_system<double, 2, 1, 1>(1e-10);
    auto Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    // Asserted at the solver, where the verdict originates ...
    auto const solved = ctrlpp::dare<double, 2, 1>(sys.A, sys.B, Q, R);
    REQUIRE_FALSE(solved.has_value());
    REQUIRE(solved.error() == ctrlpp::dare_error::non_stabilizable);

    // ... and at the gain, which forwards it rather than flattening it. Without
    // this second assertion the forwarding is untested at the surface callers
    // actually use.
    auto const result = ctrlpp::lqr_gain<double, 2, 1>(sys.A, sys.B, Q, R);
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::dare_error::non_stabilizable);
}

TEST_CASE("LQR scalar integrator analytical gain", "[lqr][hardening][precision]")
{
    Eigen::Matrix<double, 1, 1> A, B, Q, R;
    A(0, 0) = 1.0;
    B(0, 0) = 1.0;
    Q(0, 0) = 1.0;
    R(0, 0) = 1.0;

    auto result = ctrlpp::lqr_gain<double, 1, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    constexpr double eps = std::numeric_limits<double>::epsilon();

    // The analytical target is right -- the scalar Riccati fixed point for
    // A = B = Q = R = 1 is the golden ratio -- and the budget is now carried
    // there rather than guessed. The residual bounds how far the solver's answer
    // can sit from a true solution; the residual's own derivative at that
    // solution converts that into a bound on the solution itself.
    //
    // For this data the residual is r(P) = 1 - P^2 / (1 + P), so
    // r'(P) = -(P^2 + 2P) / (1 + P)^2, which is 0.854 at the golden ratio. A
    // residual bounded by the counted chain therefore bounds the forward error
    // by that chain divided by 0.854, and the gain K = P / (1 + P) contracts it
    // further by K'(P) = 1 / (1 + P)^2. Two further roundings enter on the test
    // side, one for the square root and one for the sum forming the ratio.
    double const golden = (1.0 + std::sqrt(5.0)) / 2.0;
    double const expected_K = golden / (1.0 + golden);

    double const residual_slope = (golden * golden + 2.0 * golden)
                                  / ((1.0 + golden) * (1.0 + golden));
    constexpr int analytic_ops = 2;
    double const solution_budget =
        ctrlpp::test::riccati_residual_ops<1, 1> * eps * golden / residual_slope
        + analytic_ops * eps * golden;
    double const gain_budget = solution_budget / ((1.0 + golden) * (1.0 + golden))
                               + analytic_ops * eps * expected_K;

    auto const solved = ctrlpp::dare<double, 1, 1>(A, B, Q, R);
    REQUIRE(solved.has_value());
    CAPTURE(solved->P(0, 0), golden, solution_budget);
    REQUIRE(std::abs(solved->P(0, 0) - golden) <= solution_budget);

    CAPTURE((*result)(0, 0), expected_K, gain_budget);
    REQUIRE(std::abs((*result)(0, 0) - expected_K) <= gain_budget);
}

TEST_CASE("lqr_tv_gains refuses a state sequence one shorter than the rest", "[lqr][hardening][negative]")
{
    auto const result = tv_gains_of(double_integrator_sequences(49, 50, 50, 50));
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::lqr_error::sequence_length_mismatch);
}

TEST_CASE("lqr_tv_gains refuses any one sequence one short of the others", "[lqr][hardening][negative]")
{
    constexpr std::size_t n = 50;
    for(std::size_t shortened = 0; shortened < 4; ++shortened)
    {
        auto const size_of = [&](std::size_t which) { return which == shortened ? n - 1 : n; };
        auto const result = tv_gains_of(double_integrator_sequences(size_of(0), size_of(1), size_of(2), size_of(3)));
        CAPTURE(shortened);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == ctrlpp::lqr_error::sequence_length_mismatch);
    }
}

TEST_CASE("lqr_tv_gains refuses a NaN state weighting at one step", "[lqr][hardening][negative]")
{
    auto sequences = double_integrator_sequences(20, 20, 20, 20);
    sequences.Qs[7](1, 0) = std::numeric_limits<double>::quiet_NaN();

    auto const result = tv_gains_of(sequences);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::lqr_error::non_finite_input);
}

TEST_CASE("lqr_tv_gains refuses an infinite terminal weighting", "[lqr][hardening][negative]")
{
    auto sequences = double_integrator_sequences(20, 20, 20, 20);
    sequences.Qf(0, 0) = std::numeric_limits<double>::infinity();

    auto const result = tv_gains_of(sequences);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::lqr_error::non_finite_input);
}

TEST_CASE("lqr_tv_gains on constant sequences reproduces lqr_finite bitwise", "[lqr][hardening][precision]")
{
    constexpr std::size_t horizon = 20;
    auto const sequences = double_integrator_sequences(horizon, horizon, horizon, horizon);

    // Both recursions evaluate the same expressions in the same order, so with
    // constant matrices nothing but bit-identical gains is correct.
    auto const finite = ctrlpp::lqr_finite<double, 2, 1>(sequences.As[0], sequences.Bs[0], sequences.Qs[0],
                                                         sequences.Rs[0], sequences.Qf, horizon);
    auto const result = tv_gains_of(sequences);
    REQUIRE(result.has_value());
    REQUIRE(result->size() == horizon);
    for(std::size_t k = 0; k < horizon; ++k)
    {
        CAPTURE(k);
        CHECK((*result)[k] == finite[k]);
    }
}

TEST_CASE("lqr_tv_gains on empty sequences returns no gains", "[lqr][hardening][precision]")
{
    auto const result = tv_gains_of(double_integrator_sequences(0, 0, 0, 0));
    REQUIRE(result.has_value());
    CHECK(result->empty());
}

TEST_CASE("lqr_time_varying refuses a step at the horizon", "[lqr][hardening][negative]")
{
    ctrlpp::lqr_time_varying<double, 2, 1> law(ramp_gains(50));
    input_matrix x;
    x << 1.0, -0.5;

    auto const u = law.compute(x, law.horizon());
    REQUIRE_FALSE(u.has_value());
    CHECK(u.error() == ctrlpp::lqr_error::step_out_of_range);

    auto const K = law.gain(law.horizon());
    REQUIRE_FALSE(K.has_value());
    CHECK(K.error() == ctrlpp::lqr_error::step_out_of_range);
}

TEST_CASE("lqr_time_varying serves every step before the horizon", "[lqr][hardening][precision]")
{
    auto const gains = ramp_gains(8);
    ctrlpp::lqr_time_varying<double, 2, 1> law(gains);
    input_matrix x;
    x << 1.0, -0.5;

    for(std::size_t k = 0; k < law.horizon(); ++k)
    {
        CAPTURE(k);
        auto const u = law.compute(x, k);
        auto const K = law.gain(k);
        REQUIRE(u.has_value());
        REQUIRE(K.has_value());
        CHECK(*u == (-gains[k] * x).eval());
        CHECK(*K == gains[k]);
    }
}

TEST_CASE("lqr_cost refuses two states against three inputs", "[lqr][hardening][negative]")
{
    auto const [Q, R] = trajectory_weights();
    auto const result = cost_of(dyadic_states(2), dyadic_inputs(3), Q, R);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::lqr_error::trajectory_length_mismatch);
}

TEST_CASE("lqr_cost refuses a state sequence two longer than the inputs", "[lqr][hardening][negative]")
{
    auto const [Q, R] = trajectory_weights();
    auto const result = cost_of(dyadic_states(5), dyadic_inputs(3), Q, R);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::lqr_error::trajectory_length_mismatch);
}

TEST_CASE("lqr_cost refuses a non-finite state, input or weight", "[lqr][hardening][negative]")
{
    auto const [Q, R] = trajectory_weights();
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();

    auto xs = dyadic_states(4);
    xs[2](1) = nan;
    auto us = dyadic_inputs(3);
    us[1](0) = nan;
    state_matrix Q_nan = Q;
    Q_nan(0, 1) = nan;
    weight_matrix R_nan = R;
    R_nan(0, 0) = nan;

    auto const refusals = {cost_of(xs, dyadic_inputs(3), Q, R), cost_of(dyadic_states(4), us, Q, R),
                           cost_of(dyadic_states(4), dyadic_inputs(3), Q_nan, R),
                           cost_of(dyadic_states(4), dyadic_inputs(3), Q, R_nan)};
    for(auto const& result : refusals)
    {
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == ctrlpp::lqr_error::non_finite_input);
    }
}

TEST_CASE("lqr_cost sums both accepted trajectory shapes exactly", "[lqr][hardening][precision]")
{
    auto const [Q, R] = trajectory_weights();
    auto const xs = dyadic_states(4);
    auto const us = dyadic_inputs(3);

    double stage = 0.0;
    for(std::size_t k = 0; k < us.size(); ++k)
        stage += xs[k].squaredNorm() + 2.0 * us[k].squaredNorm();
    double const terminal = xs[3].squaredNorm();

    auto const without_terminal = cost_of(dyadic_states(3), us, Q, R);
    auto const with_terminal = cost_of(xs, us, Q, R);
    REQUIRE(without_terminal.has_value());
    REQUIRE(with_terminal.has_value());
    CHECK(*without_terminal == stage);
    CHECK(*with_terminal == stage + terminal);
}
