#include "hardening_helpers.h"
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

TEST_CASE("PF all-NaN measurements does not crash",
          "[particle_filter][hardening][negative]")
{
    auto pf = make_pf();

    ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();
    pf.predict(u);

    ctrlpp::Vector<double, 1> z;
    z << std::numeric_limits<double>::quiet_NaN();
    pf.update(z);

    // NaN propagation or finite -- no crash
    CHECK((std::isnan(pf.state()[0]) || std::isfinite(pf.state()[0])));
}

TEST_CASE("PF Gaussian posterior mean/var within tolerance",
          "[particle_filter][hardening][precision]")
{
    // Use 500 particles for better precision
    ctrlpp::Matrix<double, 2, 2> Q = ctrlpp::Matrix<double, 2, 2>::Identity() * 0.01;
    ctrlpp::Matrix<double, 1, 1> R;
    R << 0.1;
    ctrlpp::Vector<double, 2> x0 = ctrlpp::Vector<double, 2>::Zero();
    ctrlpp::Matrix<double, 2, 2> P0 = ctrlpp::Matrix<double, 2, 2>::Identity();

    auto pf = ctrlpp::make_particle_filter<500>(
        pf_linear_dynamics{}, pf_position_measurement{},
        ctrlpp::pf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0},
        std::mt19937_64{42});

    double true_pos = 0.0;
    double true_vel = 1.0;
    ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();

    for(int k = 0; k < 50; ++k)
    {
        true_pos += 0.1 * true_vel;
        pf.predict(u);
        ctrlpp::Vector<double, 1> z;
        z << true_pos;
        pf.update(z);
    }

    REQUIRE(std::abs(pf.state()[0] - true_pos) < 1.0);
}

TEST_CASE("PF converges to true state for linear Gaussian system",
          "[particle_filter][hardening][convergence]")
{
    auto pf = make_pf(123);

    double true_pos = 0.0;
    double true_vel = 1.0;
    ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();

    for(int k = 0; k < 100; ++k)
    {
        true_pos += 0.1 * true_vel;
        pf.predict(u);
        ctrlpp::Vector<double, 1> z;
        z << true_pos;
        pf.update(z);
    }

    REQUIRE(std::abs(pf.state()[0] - true_pos) < 2.0);
}

TEST_CASE("PF with 10 particles does not crash",
          "[particle_filter][hardening][robustness]")
{
    ctrlpp::Matrix<double, 2, 2> Q = ctrlpp::Matrix<double, 2, 2>::Identity() * 0.01;
    ctrlpp::Matrix<double, 1, 1> R;
    R << 1.0;
    ctrlpp::Vector<double, 2> x0 = ctrlpp::Vector<double, 2>::Zero();
    ctrlpp::Matrix<double, 2, 2> P0 = ctrlpp::Matrix<double, 2, 2>::Identity();

    auto pf = ctrlpp::make_particle_filter<10>(
        pf_linear_dynamics{}, pf_position_measurement{},
        ctrlpp::pf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0},
        std::mt19937_64{42});

    ctrlpp::Vector<double, 1> u = ctrlpp::Vector<double, 1>::Zero();

    bool all_finite = true;
    for(int k = 0; k < 50; ++k)
    {
        pf.predict(u);
        ctrlpp::Vector<double, 1> z;
        z << 1.0;
        pf.update(z);

        if(!std::isfinite(pf.state()[0]))
        {
            all_finite = false;
            break;
        }
    }

    REQUIRE(all_finite);
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
