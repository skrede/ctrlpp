// What the oracles in this file decide.
//
// A particle filter's answer is a Monte-Carlo estimate, so its oracles come in
// two kinds and the file keeps them apart:
//
//  * Where the answer is DETERMINISTIC, it is asserted exactly or against a
//    counted-operation budget. A refused measurement is asserted bitwise, by
//    twin filters of which only one saw it. The weight-depletion guard is the
//    other case: when no particle explains a finite measurement the filter
//    resets to uniform weights and reports it through health, so the reported
//    state is the plain particle mean, and that is an identity rather than a
//    statistical claim. Two cases below reach that guard, one per weight
//    representation, and both assert the same identity with the same budget,
//    counted over the particle count.
//  * Where the answer is STATISTICAL, it is asserted against the exact posterior
//    the model admits. Every tracking case here is linear and Gaussian, so the
//    exact posterior is a ctrlpp::kalman_filter's, and the deviation from it is
//    the Monte-Carlo sampling error of the estimator in question -- for the mean,
//    the cloud's own standard deviation over the square root of the particle
//    count; for the variance, the relative sampling error of a second moment
//    over that many samples. The multiplier is three sampling errors, stated once
//    and used unchanged in every such case rather than tuned per case. The
//    generator seed is fixed in every case, so these assertions are
//    deterministic and the multiplier is not guarding against flakiness -- it
//    states the statistical scale of the estimator itself.
//
// Both halves of the posterior are asserted, because both are now reachable:
// `covariance()` reports the WEIGHTED second moment, which is the posterior's
// uncertainty. The unweighted dispersion of `particles()` is not a substitute
// and is not used as one -- it is measured here at 2.8 times the exact posterior
// variance and 1.9 times the exact prior, because it does not see the weights at
// all.
//
// What they deliberately do not decide. Nothing here asserts that resampling
// fired on any particular step: the event is not observable from outside. What
// IS asserted in its place is the consequence that must hold if the measurements
// conditioned the cloud at all -- the spread stays below what an unconditioned
// cloud's would necessarily have grown to.

#include "hardening_helpers.h"

#include "ctrlpp/estimation/kalman.h"
#include "ctrlpp/estimation/particle_filter.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <numbers>
#include <algorithm>

using Catch::Matchers::WithinAbs;

namespace {

// Sampling errors of slack allowed between a particle-filter estimate and the
// exact posterior. Three is the conventional bound on a Monte-Carlo estimate
// and is used unchanged in every statistical case in this file; the observed
// margins differ by an order of magnitude between cases, so a per-case value
// would be a fit rather than a statement.
constexpr double pf_standard_errors = 3.0;

/// @brief Relative sampling error of a second moment estimated from `count`
/// samples of a Gaussian: sqrt(2 / (count - 1)).
///
/// Deliberately optimistic, and the bound holds anyway. Resampling duplicates
/// particles, so the effective sample size is below the particle count and the
/// true sampling error is larger than this. Using the independent-sample formula
/// therefore makes the variance assertions HARDER to satisfy, not easier, which
/// is the direction an honest simplification has to err in.
auto second_moment_relative_sampling_error(double count) -> double
{
    return std::sqrt(2.0 / (count - 1.0));
}

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

    return ctrlpp::test::constructed(ctrlpp::make_particle_filter<100>(
        pf_linear_dynamics{}, pf_position_measurement{},
        ctrlpp::pf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0},
        std::mt19937_64{seed}));
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

struct nan_dynamics
{
    auto operator()(const ctrlpp::Vector<double, 2>&, const ctrlpp::Vector<double, 1>&) const -> ctrlpp::Vector<double, 2>
    {
        return ctrlpp::test::nan_vector<double, 2>();
    }
};

// A Gaussian policy that scores a negative measurement NaN and one above 100
// +inf, so a finite measurement can poison the likelihood on demand.
struct poisoned_likelihood
{
    double operator()(const ctrlpp::Vector<double, 1>& z, const ctrlpp::Vector<double, 1>& z_pred,
                      const ctrlpp::Matrix<double, 1, 1>& R_inv, double log_det_2piR) const
    {
        if(z(0) < 0.0)
            return std::numeric_limits<double>::quiet_NaN();
        if(z(0) > 100.0)
            return std::numeric_limits<double>::infinity();
        return ctrlpp::detail::gaussian_likelihood<double, 1>{}(z, z_pred, R_inv, log_det_2piR);
    }
};

template <typename Particles>
auto same_particles(const Particles& a, const Particles& b) -> bool
{
    for(std::size_t i = 0; i < a.size(); ++i)
    {
        if(!(a[i].array() == b[i].array()).all())
            return false;
    }
    return true;
}

template <typename Filter>
auto same_filter(const Filter& a, const Filter& b) -> bool
{
    return same_particles(a.particles(), b.particles()) && (a.state().array() == b.state().array()).all()
           && (a.covariance().array() == b.covariance().array()).all();
}

/// @brief Assert that the reported state and covariance are the plain mean and
/// second moment of the particles, as uniform weights make them.
///
/// An NP-term weighted sum accumulates at most one rounding per term, at the
/// scale of the largest operand. Uniform weights get no special case in
/// `covariance()`, so the weighted second moment reduces to the plain one with
/// the 1/NP normalization the weights carry: the measurement carried no
/// information, so the reported uncertainty is the dispersion the filter was
/// already carrying.
template <typename Filter>
void check_uniform_weight_identity(const Filter& pf)
{
    const auto particle_count = static_cast<double>(pf.particles().size());
    ctrlpp::Vector<double, 2> ref_mean = ctrlpp::Vector<double, 2>::Zero();
    for(const auto& p : pf.particles())
        ref_mean += p;
    ref_mean /= particle_count;

    const double tol = std::max(1.0, ref_mean.cwiseAbs().maxCoeff()) * particle_count * std::numeric_limits<double>::epsilon();
    CHECK(std::abs(pf.state()(0) - ref_mean(0)) <= tol);
    CHECK(std::abs(pf.state()(1) - ref_mean(1)) <= tol);

    ctrlpp::Matrix<double, 2, 2> ref_moment = ctrlpp::Matrix<double, 2, 2>::Zero();
    for(const auto& p : pf.particles())
    {
        auto d = (p - ref_mean).eval();
        ref_moment += d * d.transpose();
    }
    ref_moment /= particle_count;

    const double moment_tol = ref_moment.cwiseAbs().maxCoeff() * particle_count * std::numeric_limits<double>::epsilon();
    CHECK((pf.covariance() - ref_moment).cwiseAbs().maxCoeff() <= moment_tol);
}

void check_poisoned_likelihood_refusal(ctrlpp::weight_representation weights, double poison)
{
    ctrlpp::pf_config<double, 2, 1, 1> cfg{.Q = ctrlpp::Matrix<double, 2, 2>::Identity() * 0.01};
    cfg.weights = weights;
    const auto make = [&] {
        return ctrlpp::test::constructed(ctrlpp::make_particle_filter<100>(
            pf_linear_dynamics{}, pf_position_measurement{}, cfg, std::mt19937_64{42}, poisoned_likelihood{}));
    };
    auto a = make();
    auto b = make();
    const ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();
    ctrlpp::Vector<double, 1> poisoned;
    poisoned << poison;

    for(int k = 0; k < 10; ++k)
    {
        ctrlpp::Vector<double, 1> z;
        z << 0.1 * k;
        a.predict(u);
        b.predict(u);
        REQUIRE(a.update(z).has_value());
        REQUIRE(b.update(z).has_value());
        const auto before = a.particles();
        const auto refused = a.update(poisoned);
        REQUIRE_FALSE(refused.has_value());
        CHECK(refused.error() == ctrlpp::pf_update_error::non_finite_result);
        REQUIRE(same_particles(a.particles(), before));
    }

    REQUIRE(same_filter(a, b));
    CHECK(a.health() == ctrlpp::pf_health::ok);
}

}

TEST_CASE("PF refuses a non-finite measurement and leaves the filter bitwise unchanged",
          "[particle_filter][hardening][negative]")
{
    // Twin filters from the same seed: A sees one extra NaN measurement that B
    // never sees. If the refused step touched a weight, a particle or the
    // generator, the later resampling-capable steps would carry the difference
    // into every particle, the state and the covariance.
    auto a = make_pf();
    auto b = make_pf();
    const ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();
    const auto nan_z = ctrlpp::test::nan_vector<double, 1>();

    for(int k = 0; k < 10; ++k)
    {
        ctrlpp::Vector<double, 1> z;
        z << 0.1 * k;
        a.predict(u);
        b.predict(u);
        REQUIRE(a.update(z).has_value());
        REQUIRE(b.update(z).has_value());
        if(k == 0)
        {
            const auto refused = a.update(nan_z);
            REQUIRE_FALSE(refused.has_value());
            CHECK(refused.error() == ctrlpp::pf_update_error::non_finite_measurement);
        }
    }

    REQUIRE(same_filter(a, b));
    CHECK(a.health() == ctrlpp::pf_health::ok);
    CHECK(b.health() == ctrlpp::pf_health::ok);
}

TEST_CASE("PF reports particles a dynamics model made non-finite ahead of the measurement",
          "[particle_filter][hardening][negative]")
{
    auto pf = ctrlpp::test::constructed(ctrlpp::make_particle_filter<16>(
        nan_dynamics{}, pf_position_measurement{}, ctrlpp::pf_config<double, 2, 1, 1>{}));
    pf.predict(ctrlpp::Vector<double, 1>::Zero());

    const auto finite = pf.update(ctrlpp::Vector<double, 1>::Zero());
    REQUIRE_FALSE(finite.has_value());
    CHECK(finite.error() == ctrlpp::pf_update_error::non_finite_state);
    CHECK(pf.health() == ctrlpp::pf_health::non_finite_estimate);

    const auto non_finite = pf.update(ctrlpp::test::nan_vector<double, 1>());
    REQUIRE_FALSE(non_finite.has_value());
    CHECK(non_finite.error() == ctrlpp::pf_update_error::non_finite_state);
    CHECK(pf.health() == ctrlpp::pf_health::non_finite_estimate);
}

TEST_CASE("PF refuses a finite measurement whose likelihood is not finite and leaves the filter bitwise unchanged",
          "[particle_filter][hardening][negative]")
{
    for(const auto weights : {ctrlpp::weight_representation::log, ctrlpp::weight_representation::linear})
    {
        for(const double poison : {-1.0, 1e3})
        {
            CAPTURE(weights == ctrlpp::weight_representation::log, poison);
            check_poisoned_likelihood_refusal(weights, poison);
        }
    }
}

TEST_CASE("PF posterior mean and variance track the exact Gaussian posterior",
          "[particle_filter][hardening][precision]")
{
    // The model is linear and the noises are Gaussian, so the exact posterior is
    // the Kalman filter's and the particle estimate is a Monte-Carlo
    // approximation OF IT -- in BOTH moments. That is the reference. The bound
    // this replaces was a distance of 1.0 from the TRUTH, which is three
    // measurement standard deviations at this noise level and held for a filter
    // that had lost the target entirely; and the variance the case was named for
    // was never asserted at all.
    ctrlpp::Matrix<double, 2, 2> Q = ctrlpp::Matrix<double, 2, 2>::Identity() * 0.01;
    ctrlpp::Matrix<double, 1, 1> R;
    R << 0.1;
    ctrlpp::Vector<double, 2> x0 = ctrlpp::Vector<double, 2>::Zero();
    ctrlpp::Matrix<double, 2, 2> P0 = ctrlpp::Matrix<double, 2, 2>::Identity();

    constexpr std::size_t NP = 500;
    auto pf = ctrlpp::test::constructed(ctrlpp::make_particle_filter<NP>(
        pf_linear_dynamics{}, pf_position_measurement{},
        ctrlpp::pf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0},
        std::mt19937_64{42}));

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
        REQUIRE(pf.update(z).has_value());
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

    // The variance the case is named for, against the same exact posterior. This
    // reads the WEIGHTED second moment; the unweighted dispersion of the cloud
    // is 2.8 times the exact posterior variance here and would fail by a wide
    // margin, which is precisely why it is not a substitute for it.
    const double exact_variance = exact.covariance()(0, 0);
    const double variance_slack =
        pf_standard_errors
        * second_moment_relative_sampling_error(static_cast<double>(NP)) * exact_variance;
    CAPTURE(pf.covariance()(0, 0), exact_variance, variance_slack);
    REQUIRE(std::abs(pf.covariance()(0, 0) - exact_variance) <= variance_slack);

    // Symmetry of the reported covariance is an equality, not a tolerance: the
    // accessor symmetrizes, exactly as the other filters' covariance paths do.
    REQUIRE(pf.covariance()(0, 1) == pf.covariance()(1, 0));
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
        REQUIRE(pf.update(z).has_value());
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

    auto pf = ctrlpp::test::constructed(ctrlpp::make_particle_filter<NP>(
        pf_linear_dynamics{}, pf_position_measurement{},
        ctrlpp::pf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0},
        std::mt19937_64{42}));

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
        REQUIRE(pf.update(z).has_value());
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

TEST_CASE("PF log-mode depletion resets the weights and reports it through health",
          "[particle_filter][hardening][robustness]")
{
    auto pf = make_pf();
    pf.predict(ctrlpp::Vector<double, 1>::Zero());

    // Every squared Mahalanobis distance overflows, so every log-likelihood is
    // -inf: a finite measurement no particle explains.
    ctrlpp::Vector<double, 1> z;
    z << 1e200;
    REQUIRE(pf.update(z).has_value());
    CHECK(pf.health() == ctrlpp::pf_health::weights_reset);
    check_uniform_weight_identity(pf);
    CHECK(pf.covariance()(0, 1) == pf.covariance()(1, 0));

    const auto refused = pf.update(ctrlpp::test::nan_vector<double, 1>());
    REQUIRE_FALSE(refused.has_value());
    CHECK(refused.error() == ctrlpp::pf_update_error::non_finite_measurement);
    CHECK(pf.health() == ctrlpp::pf_health::weights_reset);
}

TEST_CASE("PF linear-mode total underflow resets the weights and reports it through health",
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

    auto pf = ctrlpp::test::constructed(ctrlpp::make_particle_filter<NP>(
        pf_linear_dynamics{}, pf_position_measurement{}, cfg, std::mt19937_64{42}));

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
    REQUIRE(pf.update(z).has_value());
    CHECK(pf.health() == ctrlpp::pf_health::weights_reset);

    // No resampling occurs (uniform weights give ESS = NP > NP/2), so the
    // extracted state must equal the plain mean of the particles.
    check_uniform_weight_identity(pf);
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

    auto pf = ctrlpp::test::constructed(ctrlpp::make_particle_filter<500>(
        bearing_dynamics{}, bearing_measurement{},
        ctrlpp::pf_config<double, 1, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0},
        std::mt19937_64{7}, bearing_wrap_likelihood{}));

    std::mt19937_64 sensor_rng{99};
    std::normal_distribution<double> noise(0.0, std::sqrt(R_scalar));
    ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();

    for(int k = 0; k < 200; ++k)
    {
        pf.predict(u);
        ctrlpp::Vector<double, 1> z;
        z << std::remainder(true_bearing + noise(sensor_rng), 2.0 * pi);
        REQUIRE(pf.update(z).has_value());
    }

    // The fused estimate should track the true bearing to within one raw
    // measurement standard deviation sqrt(R), measured as a wrapped angular error.
    double wrapped_err = std::remainder(pf.state()(0) - true_bearing, 2.0 * pi);
    CAPTURE(pf.state()(0), wrapped_err);
    CHECK(std::abs(wrapped_err) < std::sqrt(R_scalar));
}

namespace {

struct identity_dynamics
{
    auto operator()(const ctrlpp::Vector<double, 2>& x, const ctrlpp::Vector<double, 1>&) const -> ctrlpp::Vector<double, 2> { return x; }
};

using pf_config_2 = ctrlpp::pf_config<double, 2, 1, 1>;

auto refusal_of(const pf_config_2& config) -> ctrlpp::filter_error
{
    auto pf = ctrlpp::make_particle_filter<16>(identity_dynamics{}, pf_position_measurement{}, config);
    REQUIRE_FALSE(pf.has_value());
    return pf.error();
}

auto diag2(double a, double b) -> ctrlpp::Matrix<double, 2, 2>
{
    return ctrlpp::Vector<double, 2>(a, b).asDiagonal();
}

}

TEST_CASE("PF samples a rank-one process noise as the covariance it was given",
          "[particle_filter][hardening][sampling]")
{
    // The discretized white-noise-acceleration Q is rank one with an exactly zero
    // second Cholesky pivot. A point-mass prior and identity dynamics make the
    // cloud after one predict a pure sample of N(0, Q), whose entrywise
    // Monte-Carlo standard error is sqrt((Q_ij^2 + Q_ii Q_jj) / NP).
    constexpr std::size_t NP = 4000;
    ctrlpp::Matrix<double, 2, 2> Q;
    Q << 1.0, 2.0, 2.0, 4.0;
    const pf_config_2 config{.Q = Q, .R = ctrlpp::Matrix<double, 1, 1>::Identity(), .x0 = ctrlpp::Vector<double, 2>::Zero(),
                             .P0 = ctrlpp::Matrix<double, 2, 2>::Zero()};

    auto pf = ctrlpp::test::constructed(ctrlpp::make_particle_filter<NP>(identity_dynamics{}, pf_position_measurement{}, config, std::mt19937_64{42}));
    pf.predict(ctrlpp::Vector<double, 1>::Zero());

    const auto sampled = pf.covariance();
    for(int i = 0; i < 2; ++i)
    {
        for(int j = 0; j < 2; ++j)
        {
            const double standard_error = std::sqrt((Q(i, j) * Q(i, j) + Q(i, i) * Q(j, j)) / static_cast<double>(NP));
            CAPTURE(i, j, sampled(i, j), Q(i, j), standard_error);
            CHECK(std::abs(sampled(i, j) - Q(i, j)) <= pf_standard_errors * standard_error);
        }
    }
}

TEST_CASE("PF factory refuses a configuration with no Gaussian to sample or score",
          "[particle_filter][hardening][negative]")
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const pf_config_2 valid{};

    SECTION("process noise")
    {
        auto config = valid;
        config.Q = diag2(1.0, -1.0);
        CHECK(refusal_of(config) == ctrlpp::filter_error::indefinite_process_noise);
        config.Q(0, 0) = nan;
        CHECK(refusal_of(config) == ctrlpp::filter_error::non_finite_process_noise);
    }

    SECTION("measurement noise")
    {
        auto config = valid;
        config.R << -1.0;
        CHECK(refusal_of(config) == ctrlpp::filter_error::indefinite_measurement_noise);
        config.R << 0.0;
        CHECK(refusal_of(config) == ctrlpp::filter_error::singular_measurement_noise);
        config.R << nan;
        CHECK(refusal_of(config) == ctrlpp::filter_error::non_finite_measurement_noise);
    }

    SECTION("initial state and covariance")
    {
        auto config = valid;
        config.P0 = diag2(1.0, -1.0);
        CHECK(refusal_of(config) == ctrlpp::filter_error::indefinite_initial_covariance);
        config.P0(1, 1) = nan;
        CHECK(refusal_of(config) == ctrlpp::filter_error::non_finite_initial_covariance);
        config.x0(0) = nan;
        CHECK(refusal_of(config) == ctrlpp::filter_error::non_finite_initial_state);
    }

    SECTION("the earlier field in declaration order is reported")
    {
        auto config = valid;
        config.P0 = diag2(1.0, -1.0);
        config.R << 0.0;
        CHECK(refusal_of(config) == ctrlpp::filter_error::singular_measurement_noise);
        config.Q = diag2(-1.0, 1.0);
        CHECK(refusal_of(config) == ctrlpp::filter_error::indefinite_process_noise);
    }
}

TEST_CASE("PF factory refuses a singular multi-output measurement noise",
          "[particle_filter][hardening][negative]")
{
    struct full_measurement
    {
        auto operator()(const ctrlpp::Vector<double, 2>& x) const -> ctrlpp::Vector<double, 2> { return x; }
    };
    ctrlpp::pf_config<double, 2, 1, 2> config{};
    config.R = diag2(1.0, 0.0);
    auto pf = ctrlpp::make_particle_filter<16>(identity_dynamics{}, full_measurement{}, config);
    REQUIRE_FALSE(pf.has_value());
    CHECK(pf.error() == ctrlpp::filter_error::singular_measurement_noise);
}

TEST_CASE("PF factory accepts a singular positive semidefinite Q and P0",
          "[particle_filter][hardening][construction]")
{
    ctrlpp::Matrix<double, 2, 2> rank_one;
    rank_one << 1.0, 2.0, 2.0, 4.0;
    const pf_config_2 config{.Q = rank_one, .P0 = ctrlpp::Matrix<double, 2, 2>::Zero()};

    auto pf = ctrlpp::make_particle_filter<16>(identity_dynamics{}, pf_position_measurement{}, config);
    REQUIRE(pf.has_value());
    for(const auto& p : pf->particles())
        CHECK(p == config.x0);
}

TEST_CASE("PF factory samples a positive-definite configuration through its Cholesky factors",
          "[particle_filter][hardening][construction]")
{
    // Replays the draws the filter makes (one normal distribution across the
    // prior, then a fresh one per particle in each predict) through the
    // unpivoted lower Cholesky factors, so the comparison is exact equality.
    constexpr std::size_t NP = 8;
    ctrlpp::Matrix<double, 2, 2> Q;
    Q << 0.01, 0.004, 0.004, 0.02;
    ctrlpp::Matrix<double, 2, 2> P0;
    P0 << 1.0, 0.3, 0.3, 0.5;
    const pf_config_2 config{.Q = Q, .x0 = ctrlpp::Vector<double, 2>(0.3, -0.2), .P0 = P0};

    auto pf = ctrlpp::test::constructed(ctrlpp::make_particle_filter<NP>(pf_linear_dynamics{}, pf_position_measurement{}, config, std::mt19937_64{42}));
    pf.predict(ctrlpp::Vector<double, 1>::Zero());

    const ctrlpp::Matrix<double, 2, 2> L_P0 = P0.llt().matrixL();
    const ctrlpp::Matrix<double, 2, 2> L_Q = Q.llt().matrixL();
    std::mt19937_64 rng{42};
    std::array<ctrlpp::Vector<double, 2>, NP> expected{};
    std::normal_distribution<double> prior(0.0, 1.0);
    for(auto& x : expected)
    {
        ctrlpp::Vector<double, 2> noise;
        noise(0) = prior(rng);
        noise(1) = prior(rng);
        x = config.x0 + L_P0 * noise;
    }
    for(auto& x : expected)
    {
        std::normal_distribution<double> process(0.0, 1.0);
        ctrlpp::Vector<double, 2> noise;
        noise(0) = process(rng);
        noise(1) = process(rng);
        x = pf_linear_dynamics{}(x, ctrlpp::Vector<double, 1>::Zero());
        x += (L_Q * noise).eval();
    }

    for(std::size_t i = 0; i < NP; ++i)
        CHECK(pf.particles()[i] == expected[i]);
}
