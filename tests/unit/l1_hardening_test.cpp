#include "hardening_helpers.h"
#include "ctrlpp/control/l1.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <utility>

using Catch::Matchers::WithinAbs;

namespace {

using Vec1 = ctrlpp::Vector<double, 1>;

auto vec1(double v) -> Vec1
{
    Vec1 r;
    r << v;
    return r;
}

// create() is the only construction path and it is fallible, so every
// valid-input site goes through it and asserts success here. The rejection
// cases below do not use this helper: they assert the specific enumerator.
auto make_siso_controller(const ctrlpp::l1_config<double, 1, 1>& cfg, double cutoff_hz, double sample_hz)
    -> ctrlpp::l1_controller<double>
{
    auto created = ctrlpp::l1_controller<double>::create(cfg, cutoff_hz, sample_hz);
    REQUIRE(created.has_value());
    return *std::move(created);
}

auto make_siso_config() -> ctrlpp::l1_config<double, 1, 1>
{
    ctrlpp::l1_config<double, 1, 1> cfg{};
    cfg.predictor_model.A << 0.9;
    cfg.predictor_model.B << 0.5;
    cfg.predictor_model.C << 1.0;
    cfg.predictor_model.D << 0.0;
    cfg.gamma << 10.0;
    cfg.theta_min << -10.0;
    cfg.theta_max << 10.0;
    return cfg;
}

}

TEST_CASE("L1 create rejects a unit-eigenvalue predictor with singular_predictor",
          "[l1][hardening][error]")
{
    auto cfg = make_siso_config();
    cfg.predictor_model.A << 1.0; // (I - A) is exactly singular
    auto result = ctrlpp::l1_controller<double>::create(cfg, 15.0, 100.0);
    REQUIRE(!result.has_value());
    CHECK(result.error() == ctrlpp::l1_error::singular_predictor);
}

TEST_CASE("L1 create rejects a zero-input predictor with singular_dc_gain",
          "[l1][hardening][error]")
{
    auto cfg = make_siso_config();
    cfg.predictor_model.B << 0.0; // DC gain (I - A)^{-1} B is exactly zero
    auto result = ctrlpp::l1_controller<double>::create(cfg, 15.0, 100.0);
    REQUIRE(!result.has_value());
    CHECK(result.error() == ctrlpp::l1_error::singular_dc_gain);
}

TEST_CASE("L1 create rejects an overflowing feedforward gain with non_finite_gain",
          "[l1][hardening][error]")
{
    auto cfg = make_siso_config();
    // A subnormal control effectiveness makes the DC gain finite and nonzero
    // (so the singularity checks pass) while its reciprocal K_r = 1 / dc_gain
    // overflows to infinity, exercising the non-finite gain rejection.
    cfg.predictor_model.B << 1.0e-320;
    auto result = ctrlpp::l1_controller<double>::create(cfg, 15.0, 100.0);
    REQUIRE(!result.has_value());
    CHECK(result.error() == ctrlpp::l1_error::non_finite_gain);
}

TEST_CASE("L1 create rejects a filter design at the Nyquist frequency with invalid_filter_config",
          "[l1][hardening][error]")
{
    auto cfg = make_siso_config();
    auto result = ctrlpp::l1_controller<double>::create(cfg, 50.0, 100.0);
    REQUIRE(!result.has_value());
    CHECK(result.error() == ctrlpp::l1_error::invalid_filter_config);
}

TEST_CASE("L1 NaN state is rejected without touching the adaptation",
          "[l1][hardening][negative]")
{
    auto cfg = make_siso_config();
    auto ctrl = make_siso_controller(cfg, 15.0, 100.0);
    // Adapt once first, so the uncertainty estimate and the predictor state
    // hold something other than their initial values.
    REQUIRE(ctrl.evaluate(vec1(0.5), vec1(1.0)).has_value());

    const auto sigma_before = ctrl.sigma_hat();
    const auto x_hat_before = ctrl.x_hat();
    const auto x_tilde_before = ctrl.tracking_error();

    auto rejected = ctrl.evaluate(ctrlpp::test::nan_vector<double, 1>(), vec1(1.0));

    REQUIRE_FALSE(rejected.has_value());
    CHECK(rejected.error() == ctrlpp::l1_step_error::non_finite_state);
    // Exact, never a tolerance: a rejected cycle performs no arithmetic on the
    // carried state. Every piece is asserted, not a representative one.
    CHECK(sigma_before == ctrl.sigma_hat());
    CHECK(x_hat_before == ctrl.x_hat());
    CHECK(x_tilde_before == ctrl.tracking_error());
    CHECK(ctrl.health() == ctrlpp::l1_health::ok);

    // The guard is what makes the projection question moot on this path. Left
    // unguarded, the NaN would reach the adaptation, and the projection would
    // NOT sanitize it: Eigen's cwiseMax/cwiseMin return their left operand when
    // the comparison is false, and every comparison against a NaN is false. An
    // infinity is a different story -- against the finite bounds configured
    // here it would be replaced by one of them, producing a finite in-range
    // command from a meaningless estimate, which is what health() reports.
    ctrlpp::l1_controller<double> reference = make_siso_controller(cfg, 15.0, 100.0);
    REQUIRE(reference.evaluate(vec1(0.5), vec1(1.0)).has_value());

    auto after = ctrlpp::test::commanded(ctrl.evaluate(vec1(0.25), vec1(1.0)));
    auto expected = ctrlpp::test::commanded(reference.evaluate(vec1(0.25), vec1(1.0)));
    CHECK(after == expected);
    CHECK(ctrl.sigma_hat() == reference.sigma_hat());
    CHECK(ctrl.x_hat() == reference.x_hat());
}

TEST_CASE("L1 NaN reference is rejected without touching the adaptation",
          "[l1][hardening][negative]")
{
    auto cfg = make_siso_config();
    auto ctrl = make_siso_controller(cfg, 15.0, 100.0);
    REQUIRE(ctrl.evaluate(vec1(0.5), vec1(1.0)).has_value());

    const auto sigma_before = ctrl.sigma_hat();
    const auto x_hat_before = ctrl.x_hat();
    const auto x_tilde_before = ctrl.tracking_error();

    auto rejected = ctrl.evaluate(vec1(0.0), ctrlpp::test::nan_vector<double, 1>());

    REQUIRE_FALSE(rejected.has_value());
    // A bad reference names the command generator, a bad state names the sensor
    // or estimator. Different subsystems, so different enumerators.
    CHECK(rejected.error() == ctrlpp::l1_step_error::non_finite_reference);
    CHECK(sigma_before == ctrl.sigma_hat());
    CHECK(x_hat_before == ctrl.x_hat());
    CHECK(x_tilde_before == ctrl.tracking_error());
    CHECK(ctrl.health() == ctrlpp::l1_health::ok);

    ctrlpp::l1_controller<double> reference = make_siso_controller(cfg, 15.0, 100.0);
    REQUIRE(reference.evaluate(vec1(0.5), vec1(1.0)).has_value());

    auto after = ctrlpp::test::commanded(ctrl.evaluate(vec1(0.25), vec1(1.0)));
    auto expected = ctrlpp::test::commanded(reference.evaluate(vec1(0.25), vec1(1.0)));
    CHECK(after == expected);
    CHECK(ctrl.sigma_hat() == reference.sigma_hat());
    CHECK(ctrl.x_hat() == reference.x_hat());
}

TEST_CASE("L1 projection substituting a bound for an overflowed adaptation is reported",
          "[l1][hardening][negative]")
{
    // Entirely finite arguments. A large adaptation gain against a large
    // prediction error overflows the raw update to an infinity, and the
    // projection then pins it to the configured bound. The command that comes
    // out is finite and inside the output range, so nothing downstream can tell
    // it apart from a command built on a meaningful estimate. The health query
    // is the only thing that can.
    auto cfg = make_siso_config();
    cfg.gamma << std::numeric_limits<double>::max();
    auto ctrl = make_siso_controller(cfg, 15.0, 100.0);

    // x_tilde = x_hat - x = +max, and the update subtracts gamma*B'*x_tilde, so
    // the raw estimate overflows to -infinity and the projection pins it to the
    // LOWER bound.
    auto u = ctrlpp::test::commanded(ctrl.evaluate(vec1(-std::numeric_limits<double>::max()), vec1(1.0)));

    CHECK(std::isfinite(u[0]));
    CHECK(std::isfinite(ctrl.sigma_hat()[0]));
    CHECK(ctrl.sigma_hat()[0] == cfg.theta_min[0]);
    CHECK(ctrl.health() == ctrlpp::l1_health::projection_clamped_non_finite);
}

TEST_CASE("L1 zero filter bandwidth produces finite output", "[l1][hardening][negative]")
{
    auto cfg = make_siso_config();
    // Very low bandwidth filter -- should still produce finite output
    auto ctrl = make_siso_controller(cfg, 0.1, 100.0);

    auto u = ctrlpp::test::commanded(ctrl.evaluate(vec1(0.0), vec1(1.0)));
    CHECK(std::isfinite(u[0]));
}

TEST_CASE("L1 known sigma_hat after one step", "[l1][hardening][precision]")
{
    auto cfg = make_siso_config();
    cfg.gamma << 1.0;
    auto ctrl = make_siso_controller(cfg, 15.0, 100.0);

    // Step 1: x=0, r=1
    // x_hat = A_m*0 + B*(0+0) = 0
    // x_tilde = x_hat - x = 0 - 0 = 0
    // sigma_hat -= gamma * B^T * x_tilde = 0
    // sigma_hat clamped to [-10, 10] => 0
    auto u = ctrlpp::test::commanded(ctrl.evaluate(vec1(0.0), vec1(1.0)));
    REQUIRE_THAT(ctrl.sigma_hat()[0], WithinAbs(0.0, 1e-15));
    CHECK(std::isfinite(u[0]));
}

TEST_CASE("L1 tracks step reference with bounded transient", "[l1][hardening][convergence]")
{
    auto cfg = make_siso_config();
    auto ctrl = make_siso_controller(cfg, 15.0, 100.0);

    double x_plant = 0.0;
    double max_overshoot = 0.0;

    for(int k = 0; k < 1000; ++k)
    {
        auto u = ctrlpp::test::commanded(ctrl.evaluate(vec1(x_plant), vec1(1.0)));
        x_plant = 0.8 * x_plant + 0.5 * u[0];

        if(x_plant > 1.0)
            max_overshoot = std::max(max_overshoot, x_plant - 1.0);
    }

    // Should converge close to reference
    REQUIRE(std::abs(x_plant - 1.0) < 0.15);
    // Overshoot should be bounded
    CHECK(max_overshoot < 1.0);
}

TEST_CASE("L1 high gamma low bandwidth bounded output", "[l1][hardening][robustness]")
{
    auto cfg = make_siso_config();
    cfg.gamma << 1e6;
    cfg.theta_min << -100.0;
    cfg.theta_max << 100.0;
    auto ctrl = make_siso_controller(cfg, 2.0, 100.0);

    double x_plant = 0.0;
    bool all_finite = true;

    for(int k = 0; k < 500; ++k)
    {
        auto u = ctrlpp::test::commanded(ctrl.evaluate(vec1(x_plant), vec1(1.0)));
        x_plant = 0.8 * x_plant + 0.5 * u[0];

        if(!std::isfinite(u[0]) || !std::isfinite(x_plant))
        {
            all_finite = false;
            break;
        }
    }

    REQUIRE(all_finite);
}

TEST_CASE("L1 projection clamps sigma_hat within bounds", "[l1][hardening][robustness]")
{
    auto cfg = make_siso_config();
    cfg.gamma << 1000.0;
    cfg.theta_min << -5.0;
    cfg.theta_max << 5.0;
    auto ctrl = make_siso_controller(cfg, 15.0, 100.0);

    double x_plant = 0.0;

    for(int k = 0; k < 200; ++k)
    {
        auto u = ctrlpp::test::commanded(ctrl.evaluate(vec1(x_plant), vec1(1.0)));
        x_plant = 0.8 * x_plant + 0.5 * u[0];
    }

    CHECK(ctrl.sigma_hat()[0] >= -5.0 - 1e-10);
    CHECK(ctrl.sigma_hat()[0] <= 5.0 + 1e-10);
}
