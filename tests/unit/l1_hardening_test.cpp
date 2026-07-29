// What the oracles in this file decide.
//
// The controller is a state predictor, an adaptation law clamped by a
// projection, and an output low-pass filter. Every claim here is a claim about
// one of those three, and each is either exact or derived from the configured
// numbers:
//
//  * The four construction rejections name their specific enumerator and are
//    the template the rest of the file follows. They are left exactly as they
//    stand.
//  * A rejected cycle leaves every piece of carried state BITWISE unchanged and
//    a later valid cycle produces exactly what an untouched instance produces.
//  * The FIRST cycle from a zero state has a closed form with no rounding in
//    it at all: the predictor state, the prediction error and the adaptation
//    increment are exact products of zero, so the raw command is the
//    feedforward gain and the filtered command is that gain times the filter's
//    leading coefficient -- one multiplication, asserted BITWISE against the
//    coefficient the design returns.
//  * The projection is an exact clamp. `cwiseMax`/`cwiseMin` return the bound
//    itself, so the estimate is asserted against its bounds with NO slack, and
//    the cases additionally assert that a bound is ATTAINED -- without that a
//    clamp case passes on a run where the clamp never engaged.
//  * Tracking is asserted against the closed loop's own fixed point, solved by
//    hand from the configured numbers. With the plant x <- 0.8 x + 0.5 u, the
//    predictor x_hat <- 0.9 x_hat + 0.5 (u + sigma), a unit-gain filter at
//    steady state and a unit reference, a stationary point needs x_hat = x,
//    which gives x = 5 (u + sigma) and x = 2.5 u, hence sigma = -0.5 u; and
//    u = K_r - sigma with K_r = 0.2 then gives u = 0.4, x = 1 and
//    sigma = -0.2. All four are asserted, not just the plant output.
//  * Boundedness is asserted against what the projection BUYS: the estimate is
//    confined to its configured range, so the raw command is confined to the
//    feedforward term plus that range, and the filtered command to the filter's
//    bounded-input gain -- the absolute sum of the impulse response of the
//    difference equation its realized coefficients define, formed
//    independently here -- times that. No round number is involved.
//
// What they deliberately do not decide. Nothing here asserts the L1 transient
// bound from the adaptation gain and the filter bandwidth. Deriving it is its
// own piece of work with its own conditioning question, and the fixed point
// above is both exact and stronger than any envelope around it would be: it
// pins the answer rather than bracketing it. The overshoot case bounds the
// excursion by the plant's own accumulation of the command's deviation from
// that fixed point, which is a property of the plant and not of the controller.

#include "hardening_helpers.h"
#include "ctrlpp/control/l1.h"
#include "ctrlpp/dsp/biquad.h"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <utility>
#include <algorithm>

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

// The coefficients the controller's own output filter is built from, obtained
// by asking the same factory with the same arguments. Nothing is reimplemented:
// the design is the library's, and only its result is read.
auto output_filter_coefficients(double cutoff_hz, double sample_hz) -> ctrlpp::biquad_coeffs<double>
{
    auto designed = ctrlpp::biquad<double>::low_pass(cutoff_hz, sample_hz);
    REQUIRE(designed.has_value());
    return designed->coefficients();
}

// The bounded-input gain of a difference equation: the absolute sum of its
// impulse response, which is the supremum of the output a unit-bounded input
// can produce. Summed over enough samples that the remaining tail is below the
// rounding of the sum -- the slowest pole here has magnitude under 0.995, so
// twenty thousand samples leave a tail below 1e-43.
auto bounded_input_gain(const ctrlpp::biquad_coeffs<double>& c) -> double
{
    constexpr int samples = 20000;
    double w1 = 0.0;
    double w2 = 0.0;
    double sum = 0.0;
    for(int k = 0; k < samples; ++k)
    {
        const double x = (k == 0) ? 1.0 : 0.0;
        const double y = c.b0 * x + w1;
        w1 = c.b1 * x - c.a1 * y + w2;
        w2 = c.b2 * x - c.a2 * y;
        sum += std::abs(y);
    }
    return sum;
}

// The reference feedforward gain the controller forms: the reciprocal of the
// predictor's direct-current gain (I - A)^{-1} B, written out for the scalar
// predictor these cases configure.
auto feedforward_gain(const ctrlpp::l1_config<double, 1, 1>& cfg) -> double
{
    return (1.0 - cfg.predictor_model.A(0, 0)) / cfg.predictor_model.B(0, 0);
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

TEST_CASE("L1 unbounded projection rejects overflow without committing it",
          "[l1][hardening][negative]")
{
    auto cfg = make_siso_config();
    cfg.gamma << std::numeric_limits<double>::max();
    cfg.theta_min << -std::numeric_limits<double>::infinity();
    cfg.theta_max << std::numeric_limits<double>::infinity();
    auto ctrl = make_siso_controller(cfg, 15.0, 100.0);

    auto const x_hat_before = ctrl.x_hat();
    auto const sigma_before = ctrl.sigma_hat();
    auto const error_before = ctrl.tracking_error();

    auto const rejected =
        ctrl.evaluate(vec1(-std::numeric_limits<double>::max()), vec1(1.0));
    REQUIRE_FALSE(rejected.has_value());
    CHECK(rejected.error() == ctrlpp::l1_step_error::non_finite_result);
    CHECK(ctrl.x_hat() == x_hat_before);
    CHECK(ctrl.sigma_hat() == sigma_before);
    CHECK(ctrl.tracking_error() == error_before);
    CHECK(ctrl.health() == ctrlpp::l1_health::ok);

    auto reference = make_siso_controller(cfg, 15.0, 100.0);
    auto const recovered =
        ctrlpp::test::commanded(ctrl.evaluate(vec1(0.0), vec1(1.0)));
    auto const expected =
        ctrlpp::test::commanded(reference.evaluate(vec1(0.0), vec1(1.0)));
    CHECK(recovered == expected);
    CHECK(ctrl.sigma_hat() == reference.sigma_hat());
}

TEST_CASE("L1 first cycle from rest is the feedforward gain through the filter's leading coefficient",
          "[l1][hardening][precision]")
{
    // Two bandwidths four orders apart, because the claim is about the design
    // the controller was handed and not about a particular design.
    for(const double cutoff_hz : {15.0, 0.1})
    {
        auto cfg = make_siso_config();
        auto ctrl = make_siso_controller(cfg, cutoff_hz, 100.0);

        // Everything the first cycle computes is an exact product of zero: the
        // predictor propagates a zero state with a zero command and a zero
        // estimate, the prediction error is zero minus the supplied zero state,
        // and the adaptation increment is the gain times a zero. The estimate
        // is therefore exactly zero -- not near zero -- and the raw command is
        // exactly the feedforward gain times the unit reference. The filter has
        // no history, so its output is that value times its leading
        // coefficient: ONE multiplication, hence bitwise equality.
        auto u = ctrlpp::test::commanded(ctrl.evaluate(vec1(0.0), vec1(1.0)));

        const double k_r = feedforward_gain(cfg);
        const double b0 = output_filter_coefficients(cutoff_hz, 100.0).b0;

        CHECK(ctrl.sigma_hat()[0] == 0.0);
        CHECK(u[0] == b0 * k_r);
        // The predictor propagated a zero state with a zero command, so it is
        // still exactly at rest and the prediction error is exactly zero.
        CHECK(ctrl.x_hat()[0] == 0.0);
        CHECK(ctrl.tracking_error()[0] == 0.0);
    }
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

    // Every term is an exact product of zero, so the estimate is exactly zero
    // and the command is exactly the feedforward gain through the filter's
    // leading coefficient. The adaptation gain does not enter either, which is
    // why this case and the one above agree despite configuring it differently.
    CHECK(ctrl.sigma_hat()[0] == 0.0);
    CHECK(u[0] == output_filter_coefficients(15.0, 100.0).b0 * feedforward_gain(cfg));
}

TEST_CASE("L1 tracks step reference with bounded transient", "[l1][hardening][convergence]")
{
    auto cfg = make_siso_config();
    auto ctrl = make_siso_controller(cfg, 15.0, 100.0);

    // The closed loop's fixed point, solved from the configured numbers. At a
    // stationary point the prediction error is zero, so the predictor and the
    // plant agree: x = 5 (u + sigma) from the predictor's direct-current gain
    // and x = 2.5 u from the plant's, hence sigma = -u/2. The command is
    // K_r r - sigma with K_r = 0.2 and r = 1, so u = 0.2 + u/2 gives u = 0.4,
    // and then x = 1 and sigma = -0.2. The filter passes a constant unchanged,
    // its coefficients summing to unit direct-current gain by construction.
    constexpr double u_fixed = 0.4;
    constexpr double x_fixed = 1.0;
    constexpr double sigma_fixed = -0.2;
    constexpr double plant_pole = 0.8;
    constexpr double plant_input_gain = 0.5;
    constexpr double reference = 1.0;

    // Eight roundings separate the realized fixed point from the exact one:
    // the plant's two products and one sum per cycle (1-3), the controller's
    // feedforward product and the subtraction forming the raw command (4, 5),
    // and the filter's three accumulations (6-8). They do not accumulate: the
    // closed loop is a contraction, so the residue at the fixed point is one
    // cycle's worth.
    constexpr int fixed_point_ops = 8;
    const double budget = fixed_point_ops * std::numeric_limits<double>::epsilon() * x_fixed;

    double x_plant = 0.0;
    double peak_command = 0.0;
    double max_overshoot = 0.0;

    for(int k = 0; k < 1000; ++k)
    {
        auto u = ctrlpp::test::commanded(ctrl.evaluate(vec1(x_plant), vec1(reference)));
        x_plant = plant_pole * x_plant + plant_input_gain * u[0];
        peak_command = std::max(peak_command, u[0]);

        if(x_plant > reference)
            max_overshoot = std::max(max_overshoot, x_plant - reference);

        // The projection is inactive throughout, which is what makes the fixed
        // point above the one the algorithm prescribes rather than the one the
        // clamp imposes. Asserted rather than assumed.
        REQUIRE(ctrl.sigma_hat()[0] > cfg.theta_min[0]);
        REQUIRE(ctrl.sigma_hat()[0] < cfg.theta_max[0]);
    }

    // All four coordinates of the fixed point, not only the plant output.
    CHECK(std::abs(x_plant - x_fixed) <= budget);
    CHECK(std::abs(ctrl.sigma_hat()[0] - sigma_fixed) <= budget);
    CHECK(std::abs(ctrl.x_hat()[0] - x_fixed) <= budget);
    CHECK(std::abs(ctrl.tracking_error()[0]) <= budget);

    // The excursion above the reference is bounded by the PLANT: writing the
    // deviation as d <- 0.8 d + 0.5 (u - u_fixed) and summing the geometric
    // series gives a supremum of 0.5 / (1 - 0.8) = 2.5 times the largest amount
    // by which the command exceeded its own fixed-point value. That is a
    // property of the plant's accumulation and needs no number chosen here.
    const double plant_bounded_input_gain = plant_input_gain / (1.0 - plant_pole);
    CHECK(max_overshoot <= plant_bounded_input_gain * std::max(0.0, peak_command - u_fixed) + budget);
}

TEST_CASE("L1 high gamma low bandwidth bounded output", "[l1][hardening][robustness]")
{
    auto cfg = make_siso_config();
    cfg.gamma << 1e6;
    cfg.theta_min << -100.0;
    cfg.theta_max << 100.0;
    constexpr double cutoff_hz = 2.0;
    auto ctrl = make_siso_controller(cfg, cutoff_hz, 100.0);

    // What the projection buys, written out. The estimate cannot leave
    // [theta_min, theta_max], so the raw command K_r r - sigma cannot leave
    // K_r r plus that range; and a linear filter cannot produce more than its
    // bounded-input gain times the supremum of its input. Every factor comes
    // from the configuration or from the realized filter coefficients.
    const double command_bound = bounded_input_gain(output_filter_coefficients(cutoff_hz, 100.0))
        * (std::abs(feedforward_gain(cfg) * 1.0) + std::max(std::abs(cfg.theta_min[0]), std::abs(cfg.theta_max[0])));

    double x_plant = 0.0;
    bool bound_attained = false;

    for(int k = 0; k < 500; ++k)
    {
        auto u = ctrlpp::test::commanded(ctrl.evaluate(vec1(x_plant), vec1(1.0)));
        x_plant = 0.8 * x_plant + 0.5 * u[0];

        REQUIRE(std::abs(u[0]) <= command_bound);
        // No slack: the projection returns the bound itself, so the estimate is
        // inside its range exactly.
        REQUIRE(ctrl.sigma_hat()[0] >= cfg.theta_min[0]);
        REQUIRE(ctrl.sigma_hat()[0] <= cfg.theta_max[0]);
        if(ctrl.sigma_hat()[0] == cfg.theta_min[0] || ctrl.sigma_hat()[0] == cfg.theta_max[0])
            bound_attained = true;
    }

    // The clamp did engage. Without this the bound above is a claim about a run
    // in which the projection never had anything to do.
    CHECK(bound_attained);
    // The raw update never overflowed, so the estimate the bound confines is a
    // meaningful one rather than a substituted bound.
    CHECK(ctrl.health() == ctrlpp::l1_health::ok);
}

TEST_CASE("L1 projection clamps sigma_hat within bounds", "[l1][hardening][robustness]")
{
    auto cfg = make_siso_config();
    cfg.gamma << 1000.0;
    cfg.theta_min << -5.0;
    cfg.theta_max << 5.0;
    auto ctrl = make_siso_controller(cfg, 15.0, 100.0);

    double x_plant = 0.0;
    int at_lower = 0;
    int at_upper = 0;

    for(int k = 0; k < 200; ++k)
    {
        auto u = ctrlpp::test::commanded(ctrl.evaluate(vec1(x_plant), vec1(1.0)));
        x_plant = 0.8 * x_plant + 0.5 * u[0];

        // The projection is an exact clamp: cwiseMax and cwiseMin return the
        // bound itself, so no slack is admissible. A padded comparison would
        // accept a projection that overshoots its own range.
        REQUIRE(ctrl.sigma_hat()[0] >= cfg.theta_min[0]);
        REQUIRE(ctrl.sigma_hat()[0] <= cfg.theta_max[0]);

        if(ctrl.sigma_hat()[0] == cfg.theta_min[0])
            ++at_lower;
        if(ctrl.sigma_hat()[0] == cfg.theta_max[0])
            ++at_upper;
    }

    // Both bounds are reached, BITWISE, which is what proves the clamp ran at
    // all. A range assertion on a run that never saturates says nothing about
    // the projection.
    CHECK(at_lower > 0);
    CHECK(at_upper > 0);
    CHECK(ctrl.health() == ctrlpp::l1_health::ok);
}
