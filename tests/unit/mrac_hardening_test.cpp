// What the oracles in this file decide.
//
// The adaptive parameters ARE the controller's memory -- nothing re-derives
// them, they only accumulate -- so every claim here is a claim about them:
//
//  * A rejected cycle leaves both parameter matrices, the reference-model state
//    and the tracking error BITWISE unchanged, and a later valid cycle produces
//    exactly what an instance that never saw the bad sample produces.
//  * With a zero adaptation gain every increment is an exact product of zero, so
//    the parameters stay bitwise at their initial values. Exact, not near.
//  * Where the adaptation law has a closed form after one cycle, that value is
//    asserted with a budget counted from the operations that produced it.
//  * Each robustification's bound is DERIVED from its own configured parameters
//    rather than compared against a round number. Writing the update as
//    theta_{k+1} = leak_k theta_k - drive_k and taking absolute values gives the
//    comparison bound b_{k+1} = |leak_k| b_k + |drive_k|, with leak fixed by
//    sigma for sigma-modification and by delta times the error norm for
//    e-modification. That bound is horizon-independent in the sense that matters:
//    it does not accumulate where the leakage term removes memory, so an
//    implementation that dropped the leakage violates it. Measured: with the
//    leakage removed the parameter norm exceeds the sigma-modification bound by
//    sixteen orders of magnitude and the e-modification bound by a factor of
//    seventeen.
//  * The dead zone's guarantee is stated as its own threshold, not as a
//    percentage: the loop settles inside the configured band, and once inside it
//    the parameters stop moving BITWISE, which is the half a tracking bound
//    alone cannot see.
//
// What they deliberately do not decide. The comparison bound is a triangle
// inequality, so it is one-sided: it catches a missing or mis-scaled leakage
// term and it does not catch a sign error in the drive. Nothing here asserts
// convergence of the parameters to their matching values either; the reference
// model and the plant are matched, but no case claims the parameters reach the
// matching gains.

#include "hardening_helpers.h"
#include "ctrlpp/control/mrac.h"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <algorithm>

namespace {

using Vec1 = ctrlpp::Vector<double, 1>;

Vec1 vec1(double v)
{
    Vec1 r;
    r << v;
    return r;
}

auto make_ref_model() -> ctrlpp::siso_discrete_state_space<double, 1>
{
    ctrlpp::siso_discrete_state_space<double, 1> ref{};
    ref.A << 0.9;
    ref.B << 0.1;
    ref.C << 1.0;
    ref.D << 0.0;
    return ref;
}

auto make_siso_config(double gamma) -> ctrlpp::mrac_controller<double>::config_type
{
    ctrlpp::mrac_controller<double>::config_type cfg{};
    cfg.reference_model = make_ref_model();
    cfg.gamma_x << gamma;
    cfg.gamma_r << gamma;
    return cfg;
}

constexpr double eps = std::numeric_limits<double>::epsilon();

// The reference model's input coupling, which is what projects the tracking
// error onto the adaptation, and the plant the robustification cases drive.
constexpr double ref_b = 0.1;
constexpr double plant_a = 0.8;
constexpr double plant_b = 0.5;

/// Rounded operations behind one comparison between a parameter and its own
/// bound.
///
/// Enumerated: the controller forms the projection (one product), scales it by
/// the sign, by the regressor and by the gain (three more), forms the error norm
/// as a square and a square root (three more), scales the parameter by the
/// leakage (one more, two for the error-dependent form) and accumulates (one
/// more): ten. The bound recomputes the same magnitudes in a different
/// association -- the leakage factor, its product with the running bound, the
/// three products of the drive and the sum: seven. Seventeen in all, every
/// operation counted whether or not it rounds.
constexpr int robustification_bound_ops = 17;

/// Drive the plant with a robustified controller and hold each parameter to the
/// comparison bound its own configured leakage produces.
///
/// `leak_of` returns the per-cycle contraction factor from the observed error
/// norm, which is what distinguishes the two robustifications: sigma-modification
/// leaks by a constant, e-modification by a factor the error itself sets.
template <typename Controller, typename LeakFn>
void assert_parameters_within_leakage_bound(Controller& ctrl, double gamma, int steps,
                                            LeakFn leak_of)
{
    double x = 0.0;
    double bound_x = 0.0;
    double bound_r = 0.0;

    for(int k = 0; k < steps; ++k)
    {
        double const x_k = x;
        auto u = ctrlpp::test::commanded(ctrl.evaluate(vec1(x_k), vec1(1.0)));
        x = plant_a * x_k + plant_b * u[0];

        double const e = std::abs(ctrl.tracking_error()(0));
        double const leak = leak_of(e);
        bound_x = leak * bound_x + gamma * ref_b * e * std::abs(x_k);
        bound_r = leak * bound_r + gamma * ref_b * e * 1.0;

        CAPTURE(k, ctrl.theta_x()(0, 0), bound_x, ctrl.theta_r()(0, 0), bound_r);
        REQUIRE(std::abs(ctrl.theta_x()(0, 0))
                <= bound_x * (1.0 + robustification_bound_ops * eps));
        REQUIRE(std::abs(ctrl.theta_r()(0, 0))
                <= bound_r * (1.0 + robustification_bound_ops * eps));
    }
}

}

TEST_CASE("MRAC NaN state is rejected without touching the adaptive parameters",
          "[mrac][hardening][negative]")
{
    auto cfg = make_siso_config(0.5);
    ctrlpp::mrac_controller<double> ctrl(cfg);
    // Adapt once first, so the parameter matrices hold something other than
    // their initial values and "unchanged" is a real claim.
    REQUIRE(ctrl.evaluate(vec1(0.5), vec1(1.0)).has_value());

    const auto theta_x_before = ctrl.theta_x();
    const auto theta_r_before = ctrl.theta_r();
    const auto x_model_before = ctrl.x_model();
    const auto tracking_error_before = ctrl.tracking_error();

    auto rejected = ctrl.evaluate(ctrlpp::test::nan_vector<double, 1>(), vec1(1.0));

    REQUIRE_FALSE(rejected.has_value());
    CHECK(rejected.error() == ctrlpp::mrac_step_error::non_finite_state);
    // Exact, never a tolerance: a rejected cycle performs no arithmetic on the
    // carried state, so bitwise equality is the contract. Every piece is
    // asserted, not a representative one.
    CHECK(theta_x_before == ctrl.theta_x());
    CHECK(theta_r_before == ctrl.theta_r());
    CHECK(x_model_before == ctrl.x_model());
    CHECK(tracking_error_before == ctrl.tracking_error());
    CHECK(ctrl.health() == ctrlpp::mrac_health::ok);

    // The parameter matrices are the controller's memory: nothing re-derives
    // them, so had the poisoned sample been admitted, no later sample would
    // have recovered them. Compare against an instance that never saw it, and
    // compare the PARAMETERS, not only the command -- a command that happens to
    // match while the parameters diverged is the failure this catches.
    ctrlpp::mrac_controller<double> reference(cfg);
    REQUIRE(reference.evaluate(vec1(0.5), vec1(1.0)).has_value());

    auto after = ctrlpp::test::commanded(ctrl.evaluate(vec1(0.25), vec1(1.0)));
    auto expected = ctrlpp::test::commanded(reference.evaluate(vec1(0.25), vec1(1.0)));
    CHECK(after == expected);
    CHECK(ctrl.theta_x() == reference.theta_x());
    CHECK(ctrl.theta_r() == reference.theta_r());
    CHECK(ctrl.x_model() == reference.x_model());
}

TEST_CASE("MRAC rejects adaptation overflow without committing it",
          "[mrac][hardening][negative]")
{
    auto cfg = make_siso_config(std::numeric_limits<double>::max());
    ctrlpp::mrac_controller<double> controller(cfg);
    auto const theta_x_before = controller.theta_x();
    auto const theta_r_before = controller.theta_r();
    auto const model_before = controller.x_model();

    auto result = controller.evaluate(vec1(4.0), vec1(1.0));
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::mrac_step_error::non_finite_result);
    CHECK(controller.theta_x() == theta_x_before);
    CHECK(controller.theta_r() == theta_r_before);
    CHECK(controller.x_model() == model_before);
    CHECK(controller.health() == ctrlpp::mrac_health::ok);
}

TEST_CASE("MRAC NaN reference is rejected without touching the adaptive parameters",
          "[mrac][hardening][negative]")
{
    auto cfg = make_siso_config(0.5);
    ctrlpp::mrac_controller<double> ctrl(cfg);
    REQUIRE(ctrl.evaluate(vec1(0.5), vec1(1.0)).has_value());

    const auto theta_x_before = ctrl.theta_x();
    const auto theta_r_before = ctrl.theta_r();
    const auto x_model_before = ctrl.x_model();
    const auto tracking_error_before = ctrl.tracking_error();

    auto rejected = ctrl.evaluate(vec1(0.0), ctrlpp::test::nan_vector<double, 1>());

    REQUIRE_FALSE(rejected.has_value());
    // A bad reference names the command generator, a bad state names the sensor
    // or estimator. Different subsystems, so different enumerators.
    CHECK(rejected.error() == ctrlpp::mrac_step_error::non_finite_reference);
    CHECK(theta_x_before == ctrl.theta_x());
    CHECK(theta_r_before == ctrl.theta_r());
    CHECK(x_model_before == ctrl.x_model());
    CHECK(tracking_error_before == ctrl.tracking_error());
    CHECK(ctrl.health() == ctrlpp::mrac_health::ok);

    ctrlpp::mrac_controller<double> reference(cfg);
    REQUIRE(reference.evaluate(vec1(0.5), vec1(1.0)).has_value());

    auto after = ctrlpp::test::commanded(ctrl.evaluate(vec1(0.25), vec1(1.0)));
    auto expected = ctrlpp::test::commanded(reference.evaluate(vec1(0.25), vec1(1.0)));
    CHECK(after == expected);
    CHECK(ctrl.theta_x() == reference.theta_x());
    CHECK(ctrl.theta_r() == reference.theta_r());
    CHECK(ctrl.x_model() == reference.x_model());
}

TEST_CASE("MRAC zero gamma produces zero adaptation", "[mrac][hardening][negative]")
{
    auto cfg = make_siso_config(0.0);
    ctrlpp::mrac_controller<double> ctrl(cfg);

    double x = 0.0;
    for(int k = 0; k < 100; ++k)
    {
        auto u = ctrlpp::test::commanded(ctrl.evaluate(vec1(x), vec1(1.0)));
        x = 0.8 * x + 0.5 * u[0];
    }

    // With a zero adaptation gain every increment is an exact product of zero, so
    // the parameters are still bitwise at their initial values after a hundred
    // cycles. A tolerance said less than the arithmetic guarantees.
    REQUIRE(ctrl.theta_x()(0, 0) == 0.0);
    REQUIRE(ctrl.theta_r()(0, 0) == 0.0);
}

TEST_CASE("MRAC huge gamma 1e15 commands zero and then the adapted gain",
          "[mrac][hardening][negative]")
{
    auto cfg = make_siso_config(1e15);
    ctrlpp::mrac_controller<double> ctrl(cfg);

    // First cycle: the state is zero and both parameters are still at their zero
    // initial values when the command is formed, so the command is exactly zero
    // whatever the gain is. Finiteness could not see the adaptation at all --
    // which is the point of the case, since an adaptation gain of 1e15 is only
    // dangerous from the SECOND cycle onwards.
    auto u = ctrlpp::test::commanded(ctrl.evaluate(vec1(0.0), vec1(1.0)));
    CHECK(u[0] == 0.0);

    // By then theta_r carries one full increment: the reference model has moved
    // to 0.1, the tracking error is -0.1, the projection through B is -0.01, and
    // the gain multiplies that by 1e15. With the state still zero the command is
    // theta_r alone.
    double const e_proj = ref_b * (0.0 - ref_b);
    double const expected = -e_proj * 1e15;
    constexpr int adaptation_ops = 3;

    auto u2 = ctrlpp::test::commanded(ctrl.evaluate(vec1(0.0), vec1(1.0)));
    CAPTURE(u2[0], expected);
    CHECK(std::abs(u2[0] - expected) <= adaptation_ops * eps * std::abs(expected));
    CHECK(ctrl.theta_x()(0, 0) == 0.0);
}

TEST_CASE("MRAC known first-order gain after one step", "[mrac][hardening][precision]")
{
    auto cfg = make_siso_config(0.5);
    ctrlpp::mrac_controller<double> ctrl(cfg);

    // Step 1: x=0, r=1
    // x_model advances: 0.9*0 + 0.1*1 = 0.1
    // tracking_error = x - x_model = 0 - 0.1 = -0.1
    // e_proj = B^T * tracking_error = 0.1 * (-0.1) = -0.01
    // theta_x -= sign_b * e_proj * x^T * gamma_x = 1 * (-0.01) * 0 * 0.5 = 0
    // theta_r -= sign_b * e_proj * r^T * gamma_r = 1 * (-0.01) * 1 * 0.5 = 0.005
    REQUIRE(ctrl.evaluate(vec1(0.0), vec1(1.0)).has_value());

    // theta_x is exact: its increment carries the state as a factor and the state
    // is exactly zero.
    REQUIRE(ctrl.theta_x()(0, 0) == 0.0);

    // theta_r carries exactly two roundings against the decimal literal. Only the
    // projection product rounds inside the controller -- the model advance is an
    // exact copy of the input coupling, the difference from zero is exact, the
    // sign and the reference are unit factors, and the adaptation gain of one
    // half scales exactly in a binary radix -- and the literal below rounds once
    // more.
    constexpr int adaptation_ops = 2;
    CAPTURE(ctrl.theta_r()(0, 0));
    REQUIRE(std::abs(ctrl.theta_r()(0, 0) - 0.005) <= adaptation_ops * eps * 0.005);
}

TEST_CASE("MRAC with a dead zone settles inside its own band and then stops adapting",
          "[mrac][hardening][convergence]")
{
    using DeadZoneMrac = ctrlpp::mrac_controller<double, 1, 1, ctrlpp::dead_zone>;
    DeadZoneMrac::config_type cfg{};
    cfg.reference_model = make_ref_model();
    cfg.gamma_x << 0.5;
    cfg.gamma_r << 0.5;
    constexpr double threshold = 0.01;
    cfg.robustification.threshold = threshold;

    DeadZoneMrac ctrl(cfg);

    constexpr int steps = 5000;
    constexpr int settled_from = 2500;
    double x = 0.0;
    double theta_x_settled = 0.0;
    double theta_r_settled = 0.0;

    for(int k = 0; k < steps; ++k)
    {
        auto u = ctrlpp::test::commanded(ctrl.evaluate(vec1(x), vec1(1.0)));
        x = plant_a * x + plant_b * u[0];

        if(k == settled_from)
        {
            theta_x_settled = ctrl.theta_x()(0, 0);
            theta_r_settled = ctrl.theta_r()(0, 0);
        }
        else if(k > settled_from)
        {
            // The dead zone's defining behavior: once the weighted error norm is
            // inside the band the adaptation is not merely small, it does not run
            // at all, so both parameters are bitwise frozen. Measured: the
            // adaptation stops permanently from cycle 62. A percentage tracking
            // bound could not see this, and it is what an implementation with the
            // dead zone removed would fail immediately.
            REQUIRE(ctrl.theta_x()(0, 0) == theta_x_settled);
            REQUIRE(ctrl.theta_r()(0, 0) == theta_r_settled);
        }
    }

    // The tracking bound is the configured threshold itself, in the same weighted
    // norm the policy tests, not a round percentage of the model state. The
    // adaptation runs exactly while the error is outside the band, so a settled
    // loop is one whose error is inside it. Measured: 7.08e-3 against a threshold
    // of 1e-2, where the five percent it replaces allowed 5e-2.
    auto const x_m = ctrl.x_model()[0];
    CAPTURE(x, x_m, std::abs(x - x_m));
    REQUIRE(std::abs(x - x_m) <= threshold);
}

TEST_CASE("MRAC sigma-modification holds the parameters inside its own leakage bound",
          "[mrac][hardening][robustness]")
{
    using SigmaMrac = ctrlpp::mrac_controller<double, 1, 1, ctrlpp::sigma_modification>;
    SigmaMrac::config_type cfg{};
    cfg.reference_model = make_ref_model();
    constexpr double gamma = 10.0;
    constexpr double sigma = 1.0;
    cfg.gamma_x << gamma;
    cfg.gamma_r << gamma;
    cfg.robustification.sigma = sigma;

    SigmaMrac ctrl(cfg);

    // Sigma-modification contracts the parameter by a constant factor every
    // cycle, so the bound's leakage is |1 - sigma| and does not depend on the
    // signals. At the configured sigma of one the contraction is total: the
    // parameter carries no memory beyond a single cycle's drive, which is why the
    // bound below is tight rather than generous. Measured: the realized parameter
    // sits at the bound to within two units in the last place, and the same bound
    // is exceeded by a factor of 1.8e16 when the leakage term is removed.
    assert_parameters_within_leakage_bound(ctrl, gamma, 2000,
                                           [](double) { return std::abs(1.0 - sigma); });
}

TEST_CASE("MRAC e-modification holds the parameters inside its own leakage bound",
          "[mrac][hardening][robustness]")
{
    using EmodMrac = ctrlpp::mrac_controller<double, 1, 1, ctrlpp::e_modification>;
    EmodMrac::config_type cfg{};
    cfg.reference_model = make_ref_model();
    constexpr double gamma = 10.0;
    constexpr double delta = 1.0;
    cfg.gamma_x << gamma;
    cfg.gamma_r << gamma;
    cfg.robustification.delta = delta;

    EmodMrac ctrl(cfg);

    // e-modification's leakage is proportional to the tracking error, so it
    // fades as the loop settles and the bound has to be carried cycle by cycle
    // rather than collapsed into a single ultimate value. Collapsing it -- taking
    // the worst leakage over the run and dividing the largest drive by one minus
    // it -- gives a bound 24000 times above the realized parameter, which is a
    // derived-looking number that could not fail. The per-cycle recursion is
    // within one unit in the last place instead, and the same bound is exceeded
    // by a factor of seventeen when the leakage term is removed.
    assert_parameters_within_leakage_bound(
        ctrl, gamma, 2000, [](double e) { return std::abs(1.0 - delta * e); });
}
