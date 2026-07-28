// What the oracles in this file decide.
//
// A particle filter's answer is a Monte-Carlo estimate, so its oracles come in
// two kinds and the file keeps them apart:
//
//  * Where the answer is DETERMINISTIC, it is asserted exactly or against a
//    counted-operation budget. The weight-degeneracy guard is the important
//    case: when every log-weight is non-finite the filter resets to uniform
//    weights, so the reported state is the plain particle mean, and that is an
//    identity rather than a statistical claim. Two cases below reach that guard
//    from opposite directions -- a non-finite measurement and a total underflow
//    -- and both assert the same identity with the same budget, counted over the
//    particle count.
//  * Where the answer is STATISTICAL, it is asserted against the exact posterior
//    the model admits. Every tracking case here is linear and Gaussian, so the
//    exact posterior is a ctrlpp::kalman_filter's, and the deviation from it is
//    the Monte-Carlo standard error of a weighted mean over the particle count:
//    the cloud's own sample standard deviation over the square root of that
//    count. The multiplier is three standard errors, stated once and used
//    unchanged in every such case rather than tuned per case. The generator seed
//    is fixed in every case, so these assertions are deterministic and the
//    multiplier is not guarding against flakiness -- it states the statistical
//    scale of the estimator itself.
//
// What they deliberately do not decide. The posterior COVARIANCE is not
// asserted anywhere, and that is a gap rather than a choice: the public surface
// exposes the particles and the weighted mean but neither the weights nor a
// posterior covariance, and the unweighted cloud variance is not a proxy for
// either -- it is measured here at nearly three times the exact posterior
// variance and about twice the prior. Nothing here asserts that resampling
// fired on any particular step either, for the same reason: the event is not
// observable from outside. What IS asserted in its place is the consequence
// that must hold if the measurements conditioned the cloud at all -- the spread
// stays below what an unconditioned cloud's would necessarily have grown to.

#include "hardening_helpers.h"

#include "ctrlpp/estimation/kalman.h"
#include "ctrlpp/estimation/particle_filter.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <random>
#include <numbers>
#include <algorithm>

using Catch::Matchers::WithinAbs;

namespace {

// Standard errors of slack allowed between a particle-filter estimate and the
// exact posterior mean. Three is the conventional bound on a Monte-Carlo mean
// and is used unchanged in every statistical case in this file; the observed
// margins differ by an order of magnitude between cases, so a per-case value
// would be a fit rather than a statement.
constexpr double pf_standard_errors = 3.0;

/// @brief Sample mean and standard deviation of the particle cloud in one
/// coordinate, and the Monte-Carlo standard error of a mean over that cloud.
///
/// The standard error is the cloud's own dispersion over the square root of the
/// particle count -- derived from the filter's realized state, not from the
/// measurement-noise magnitude, which would only coincide with it when the
/// posterior happens to be as wide as one observation.
template <typename Filter>
auto cloud_statistics(const Filter& pf, int coordinate)
{
    struct result
    {
        double standard_deviation;
        double standard_error;
    };
    const auto& particles = pf.particles();
    const auto count = static_cast<double>(particles.size());

    double mean = 0.0;
    for(const auto& p : particles)
        mean += p(coordinate);
    mean /= count;

    double variance = 0.0;
    for(const auto& p : particles)
        variance += (p(coordinate) - mean) * (p(coordinate) - mean);
    variance /= count - 1.0;

    const double deviation = std::sqrt(variance);
    return result{deviation, deviation / std::sqrt(count)};
}

/// @brief The linear system `pf_linear_dynamics` and the position measurement
/// together ARE, so the exact posterior can be computed alongside.
auto equivalent_system() -> ctrlpp::discrete_state_space<double, 2, 1, 1>
{
    ctrlpp::discrete_state_space<double, 2, 1, 1> sys;
    sys.A << 1.0, 0.1, 0.0, 1.0;
    // The particle cases all drive the filter with a zero input, so the input
    // column is irrelevant to the posterior and is left at zero rather than
    // being given a value the reference would never use.
    sys.B << 0.0, 0.0;
    sys.C << 1.0, 0.0;
    sys.D << 0.0;
    return sys;
}

struct pf_linear_dynamics
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

struct pf_position_measurement
{
    auto operator()(const ctrlpp::Vector<double, 2>& x) const -> ctrlpp::Vector<double, 1>
    {
        ctrlpp::Vector<double, 1> z;
        z(0) = x(0);
        return z;
    }
};

auto make_pf(std::size_t seed = 42)
{
    ctrlpp::Matrix<double, 2, 2> Q = ctrlpp::Matrix<double, 2, 2>::Identity() * 0.01;
    ctrlpp::Matrix<double, 1, 1> R;
    R << 0.5;
    ctrlpp::Vector<double, 2> x0 = ctrlpp::Vector<double, 2>::Zero();
    ctrlpp::Matrix<double, 2, 2> P0 = ctrlpp::Matrix<double, 2, 2>::Identity();

    return ctrlpp::make_particle_filter<100>(
        pf_linear_dynamics{}, pf_position_measurement{},
        ctrlpp::pf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0},
        std::mt19937_64{seed});
}

// Scalar random-walk bearing model: the state is an unwrapped angle advanced by
// the input rate; the sensor reports the same angle (wrapped by the sensor
// stage in the test, not the model), so the +/-pi cut is handled entirely by
// the likelihood policy.
struct bearing_dynamics
{
    auto operator()(const ctrlpp::Vector<double, 1>& x,
                    const ctrlpp::Vector<double, 1>& u) const -> ctrlpp::Vector<double, 1>
    {
        ctrlpp::Vector<double, 1> xn;
        xn(0) = x(0) + u(0);
        return xn;
    }
};

struct bearing_measurement
{
    auto operator()(const ctrlpp::Vector<double, 1>& x) const -> ctrlpp::Vector<double, 1>
    {
        return x;
    }
};

// Bearing likelihood policy: wraps the scalar innovation into [-pi, pi] via
// std::remainder before the Gaussian, so an innovation straddling the +/-pi cut
// is scored by its true angular separation rather than the naive linear
// difference. Matches the pf_likelihood_model signature.
struct bearing_wrap_likelihood
{
    double operator()(const ctrlpp::Vector<double, 1>& z, const ctrlpp::Vector<double, 1>& z_pred,
                      const ctrlpp::Matrix<double, 1, 1>& R_inv, double log_det_2piR) const
    {
        double d = std::remainder(z(0) - z_pred(0), 2.0 * std::numbers::pi_v<double>);
        double mahal = d * R_inv(0, 0) * d;
        return -0.5 * mahal - 0.5 * log_det_2piR;
    }
};

}

TEST_CASE("PF non-finite measurement recovers to the uniform particle mean",
          "[particle_filter][hardening][negative]")
{
    auto pf = make_pf();

    ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();
    pf.predict(u);

    ctrlpp::Vector<double, 1> z;
    z << std::numeric_limits<double>::quiet_NaN();
    pf.update(z);

    // The disjunction this replaces accepted both answers and so decided
    // nothing, and it hid a contract the filter genuinely has. A non-finite
    // measurement makes every log-weight non-finite; the normalization tests the
    // maximum log-weight for finiteness and, failing it, RESETS to uniform
    // weights. Uniform weights give an effective sample size equal to the
    // particle count, which is above the resampling threshold, so no resampling
    // occurs and the reported state is the plain mean of the particles -- FINITE,
    // not NaN.
    //
    // This is the identical property, reached through the identical guard, that
    // "PF linear-mode total underflow recovers to the uniform particle mean"
    // asserts further down this file; the two differ only in how the weights were
    // destroyed. Both use the same budget: an NP-term weighted sum
    // accumulates at most one rounding per term, at the scale of the largest
    // operand.
    const auto particle_count = static_cast<double>(pf.particles().size());
    ctrlpp::Vector<double, 2> ref_mean = ctrlpp::Vector<double, 2>::Zero();
    for(const auto& p : pf.particles())
        ref_mean += p;
    ref_mean /= particle_count;

    const double tol = std::max(1.0, ref_mean.cwiseAbs().maxCoeff()) * particle_count
                       * std::numeric_limits<double>::epsilon();

    CHECK(std::isfinite(pf.state()[0]));
    CHECK(std::isfinite(pf.state()[1]));
    CHECK(std::abs(pf.state()(0) - ref_mean(0)) <= tol);
    CHECK(std::abs(pf.state()(1) - ref_mean(1)) <= tol);
}

TEST_CASE("PF posterior mean tracks the exact Gaussian posterior",
          "[particle_filter][hardening][precision]")
{
    // The model is linear and the noises are Gaussian, so the exact posterior is
    // the Kalman filter's and the particle estimate is a Monte-Carlo
    // approximation OF IT. That is the reference. The bound this replaces was a
    // distance of 1.0 from the TRUTH, which is three measurement standard
    // deviations at this noise level and held for a filter that had lost the
    // target entirely.
    ctrlpp::Matrix<double, 2, 2> Q = ctrlpp::Matrix<double, 2, 2>::Identity() * 0.01;
    ctrlpp::Matrix<double, 1, 1> R;
    R << 0.1;
    ctrlpp::Vector<double, 2> x0 = ctrlpp::Vector<double, 2>::Zero();
    ctrlpp::Matrix<double, 2, 2> P0 = ctrlpp::Matrix<double, 2, 2>::Identity();

    constexpr std::size_t NP = 500;
    auto pf = ctrlpp::make_particle_filter<NP>(
        pf_linear_dynamics{}, pf_position_measurement{},
        ctrlpp::pf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0},
        std::mt19937_64{42});

    const auto sys = equivalent_system();
    auto exact = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(
        sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}));

    double true_pos = 0.0;
    double true_vel = 1.0;
    ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();
    Eigen::Matrix<double, 1, 1> u_exact = Eigen::Matrix<double, 1, 1>::Zero();

    for(int k = 0; k < 50; ++k)
    {
        true_pos += 0.1 * true_vel;
        pf.predict(u);
        exact.predict(u_exact);
        ctrlpp::Vector<double, 1> z;
        z << true_pos;
        pf.update(z);
        Eigen::Matrix<double, 1, 1> z_exact;
        z_exact << true_pos;
        REQUIRE(exact.update(z_exact).has_value());
    }

    // The generator seed is fixed, so this comparison is deterministic; the
    // multiplier states the estimator's statistical scale rather than absorbing
    // run-to-run variation.
    const auto stats = cloud_statistics(pf, 0);
    CAPTURE(pf.state()[0], exact.state()[0], stats.standard_error);
    REQUIRE(std::abs(pf.state()[0] - exact.state()[0])
            <= pf_standard_errors * stats.standard_error);
}

TEST_CASE("PF converges to the exact posterior for a linear Gaussian system",
          "[particle_filter][hardening][convergence]")
{
    auto pf = make_pf(123);

    ctrlpp::Matrix<double, 2, 2> Q = ctrlpp::Matrix<double, 2, 2>::Identity() * 0.01;
    ctrlpp::Matrix<double, 1, 1> R;
    R << 0.5;
    ctrlpp::Vector<double, 2> x0 = ctrlpp::Vector<double, 2>::Zero();
    ctrlpp::Matrix<double, 2, 2> P0 = ctrlpp::Matrix<double, 2, 2>::Identity();

    const auto sys = equivalent_system();
    auto exact = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(
        sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}));

    double true_pos = 0.0;
    double true_vel = 1.0;
    ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();
    Eigen::Matrix<double, 1, 1> u_exact = Eigen::Matrix<double, 1, 1>::Zero();

    for(int k = 0; k < 100; ++k)
    {
        true_pos += 0.1 * true_vel;
        pf.predict(u);
        exact.predict(u_exact);
        ctrlpp::Vector<double, 1> z;
        z << true_pos;
        pf.update(z);
        Eigen::Matrix<double, 1, 1> z_exact;
        z_exact << true_pos;
        REQUIRE(exact.update(z_exact).has_value());
    }

    const auto stats = cloud_statistics(pf, 0);
    CAPTURE(pf.state()[0], exact.state()[0], stats.standard_error);
    REQUIRE(std::abs(pf.state()[0] - exact.state()[0])
            <= pf_standard_errors * stats.standard_error);
}

TEST_CASE("PF with 10 particles still conditions on its measurements",
          "[particle_filter][hardening][robustness]")
{
    constexpr std::size_t NP = 10;
    constexpr int steps = 50;
    ctrlpp::Matrix<double, 2, 2> Q = ctrlpp::Matrix<double, 2, 2>::Identity() * 0.01;
    ctrlpp::Matrix<double, 1, 1> R;
    R << 1.0;
    ctrlpp::Vector<double, 2> x0 = ctrlpp::Vector<double, 2>::Zero();
    ctrlpp::Matrix<double, 2, 2> P0 = ctrlpp::Matrix<double, 2, 2>::Identity();

    auto pf = ctrlpp::make_particle_filter<NP>(
        pf_linear_dynamics{}, pf_position_measurement{},
        ctrlpp::pf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0},
        std::mt19937_64{42});

    const auto sys = equivalent_system();
    auto exact = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(
        sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}));

    ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();
    Eigen::Matrix<double, 1, 1> u_exact = Eigen::Matrix<double, 1, 1>::Zero();

    for(int k = 0; k < steps; ++k)
    {
        pf.predict(u);
        exact.predict(u_exact);
        ctrlpp::Vector<double, 1> z;
        z << 1.0;
        pf.update(z);
        Eigen::Matrix<double, 1, 1> z_exact;
        z_exact << 1.0;
        REQUIRE(exact.update(z_exact).has_value());
    }

    // Ten particles is a depletion stress case, and the Monte-Carlo standard
    // error is correspondingly large -- the estimate sits well over one standard
    // error from the exact posterior here, which is why the bound had to be
    // checked by running before it could be asserted rather than predicted.
    const auto stats = cloud_statistics(pf, 0);
    CAPTURE(pf.state()[0], exact.state()[0], stats.standard_error);
    REQUIRE(std::abs(pf.state()[0] - exact.state()[0])
            <= pf_standard_errors * stats.standard_error);

    // That the measurements conditioned the cloud at all is a separate claim, and
    // it has a falsifier that does not need the resampling event to be
    // observable: a cloud propagated without any conditioning accumulates the
    // process noise every step, so its position variance would be at least
    // P0(0,0) + steps * Q(0,0) -- and that ignores the velocity uncertainty it
    // also integrates, so it is a strict lower bound on the unconditioned spread.
    const double unconditioned_floor = std::sqrt(P0(0, 0) + steps * Q(0, 0));
    CAPTURE(stats.standard_deviation, unconditioned_floor);
    REQUIRE(stats.standard_deviation < unconditioned_floor);
}

TEST_CASE("PF linear-mode total underflow recovers to the uniform particle mean",
          "[particle_filter][hardening][robustness]")
{
    ctrlpp::Matrix<double, 2, 2> Q = ctrlpp::Matrix<double, 2, 2>::Identity() * 0.01;
    ctrlpp::Matrix<double, 1, 1> R;
    const double R_scalar = 0.5;
    R << R_scalar;
    ctrlpp::Vector<double, 2> x0 = ctrlpp::Vector<double, 2>::Zero();
    ctrlpp::Matrix<double, 2, 2> P0 = ctrlpp::Matrix<double, 2, 2>::Identity();

    constexpr std::size_t NP = 100;
    ctrlpp::pf_config<double, 2, 1, 1> cfg{.Q = Q, .R = R, .x0 = x0, .P0 = P0};
    cfg.weights = ctrlpp::weight_representation::linear;

    auto pf = ctrlpp::make_particle_filter<NP>(
        pf_linear_dynamics{}, pf_position_measurement{}, cfg, std::mt19937_64{42});

    // exp(x) underflows to zero once x < log(denorm_min); the Gaussian
    // log-likelihood carries a -0.5 * (z - x_i)^2 / R term, so any innovation
    // whose squared Mahalanobis distance exceeds -2 * log(denorm_min) underflows.
    // Place z several underflow radii beyond the furthest particle so every
    // particle's weight underflows to zero (total underflow).
    const double underflow_mahal = -2.0 * std::log(std::numeric_limits<double>::denorm_min());
    const double underflow_radius = std::sqrt(underflow_mahal * R_scalar);
    double max_pos = 0.0;
    for(const auto& p : pf.particles())
        max_pos = std::max(max_pos, std::abs(p(0)));

    ctrlpp::Vector<double, 1> z;
    z << max_pos + 2.0 * underflow_radius;
    pf.update(z);

    // No resampling occurs (uniform weights give ESS = NP > NP/2), so the
    // extracted state must be finite and equal the plain mean of the particles.
    ctrlpp::Vector<double, 2> ref_mean = ctrlpp::Vector<double, 2>::Zero();
    for(const auto& p : pf.particles())
        ref_mean += p;
    ref_mean /= static_cast<double>(NP);

    // An NP-term weighted sum accumulates at most O(NP) ulps of rounding.
    const double tol = std::max(1.0, ref_mean.cwiseAbs().maxCoeff())
                       * static_cast<double>(NP) * std::numeric_limits<double>::epsilon();

    CHECK(std::isfinite(pf.state()(0)));
    CHECK(std::isfinite(pf.state()(1)));
    CHECK(std::abs(pf.state()(0) - ref_mean(0)) <= tol);
    CHECK(std::abs(pf.state()(1) - ref_mean(1)) <= tol);
}

TEST_CASE("PF bearing wrap likelihood fixes the +/-pi cut where the Gaussian policy inverts",
          "[particle_filter][hardening][likelihood]")
{
    const double pi = std::numbers::pi_v<double>;
    const double R_scalar = 0.1;
    ctrlpp::Matrix<double, 1, 1> R_inv;
    R_inv << 1.0 / R_scalar;
    const double log_det_2piR = std::log(2.0 * pi * R_scalar);

    bearing_wrap_likelihood wrap;
    ctrlpp::detail::gaussian_likelihood<double, 1> gauss;

    // A measurement at -3.1 rad against a prediction at +3.1 rad is only ~0.083
    // rad apart across the cut, yet 6.2 rad apart as a naive linear difference.
    ctrlpp::Vector<double, 1> z_neg;
    z_neg << -3.1;
    ctrlpp::Vector<double, 1> zp_pos;
    zp_pos << 3.1;
    ctrlpp::Vector<double, 1> zp_zero;
    zp_zero << 0.0;

    double wrap_across = wrap(z_neg, zp_pos, R_inv, log_det_2piR);
    double wrap_far = wrap(z_neg, zp_zero, R_inv, log_det_2piR);
    double gauss_across = gauss(z_neg, zp_pos, R_inv, log_det_2piR);
    double gauss_far = gauss(z_neg, zp_zero, R_inv, log_det_2piR);

    CAPTURE(wrap_across, wrap_far, gauss_across, gauss_far);
    // Wrap policy: the across-cut pair is the more likely one.
    CHECK(wrap_across > wrap_far);
    // Default Gaussian policy inverts the ordering on the same inputs.
    CHECK(gauss_across < gauss_far);

    // PF-level: a stationary bearing target sitting on the +/-pi cut. Its wrapped
    // sensor readings straddle both signs; only the wrap likelihood fuses them
    // to the true angle.
    const double true_bearing = 3.1;
    ctrlpp::Matrix<double, 1, 1> Q;
    Q << 0.001;
    ctrlpp::Matrix<double, 1, 1> R;
    R << R_scalar;
    ctrlpp::Vector<double, 1> x0;
    x0 << 3.0;
    ctrlpp::Matrix<double, 1, 1> P0;
    P0 << 0.25;

    auto pf = ctrlpp::make_particle_filter<500>(
        bearing_dynamics{}, bearing_measurement{},
        ctrlpp::pf_config<double, 1, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0},
        std::mt19937_64{7}, bearing_wrap_likelihood{});

    std::mt19937_64 sensor_rng{99};
    std::normal_distribution<double> noise(0.0, std::sqrt(R_scalar));
    ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();

    for(int k = 0; k < 200; ++k)
    {
        pf.predict(u);
        ctrlpp::Vector<double, 1> z;
        z << std::remainder(true_bearing + noise(sensor_rng), 2.0 * pi);
        pf.update(z);
    }

    // The fused estimate should track the true bearing to within one raw
    // measurement standard deviation sqrt(R), measured as a wrapped angular error.
    double wrapped_err = std::remainder(pf.state()(0) - true_bearing, 2.0 * pi);
    CAPTURE(pf.state()(0), wrapped_err);
    CHECK(std::abs(wrapped_err) < std::sqrt(R_scalar));
}
