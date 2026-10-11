#include "hardening_helpers.h"
#include "ctrlpp/pid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <random>
#include <vector>
#include <cstddef>
#include <numbers>

using Catch::Matchers::WithinAbs;

namespace {

using Vec1 = ctrlpp::Vector<double, 1>;

constexpr double tol = 1e-12;

Vec1 vec1(double v) { Vec1 r; r << v; return r; }

template <typename Scalar, typename PA>
auto detector_config(Scalar amplitude, Scalar ultimate_frequency) -> typename ctrlpp::pid<Scalar, 1, PA>::config_type
{
    typename ctrlpp::pid<Scalar, 1, PA>::config_type cfg{};
    cfg.kp = ctrlpp::Vector<Scalar, 1>::Constant(Scalar{1});
    cfg.template policy<PA>().amplitude = {amplitude};
    cfg.template policy<PA>().ultimate_frequency = {ultimate_frequency};
    return cfg;
}

// Every long-run case below samples at h = 1/100 s with omega_u = 2 pi / 10, so T_u = 10 s and the
// default window is T_sup = 50 T_u = 500 s.
constexpr int steps_per_period = 1000;

template <typename Scalar>
auto ultimate_frequency() -> Scalar
{
    return Scalar{2} * std::numbers::pi_v<Scalar> / Scalar{10};
}

struct detector_run
{
    std::vector<int> sign_changes;
    int first_detection;
    bool held_after_first_detection;
    double final_metric;
};

// The error is fed as the setpoint against a zero measurement with kp = 1, so the controller's error is
// the signal itself. Sign changes are counted on that sequence, an exact zero being no sign.
template <typename Scalar, typename Signal>
auto run_detector(Scalar amplitude, int steps, Signal&& signal) -> detector_run
{
    using PA = ctrlpp::perf_assessment<ctrlpp::oscillation_detect>;
    using Pid = ctrlpp::pid<Scalar, 1, PA>;
    using vector_t = typename Pid::vector_t;
    auto pid = ctrlpp::test::constructed(Pid::create(detector_config<Scalar, PA>(amplitude, ultimate_frequency<Scalar>())));
    detector_run run{{}, -1, true, 0.0};
    Scalar previous{0};
    for(int k = 0; k < steps; ++k)
    {
        const Scalar e = signal(k);
        if(e != Scalar{0} && previous != Scalar{0} && (e > Scalar{0}) != (previous > Scalar{0}))
            run.sign_changes.push_back(k);
        if(e != Scalar{0})
            previous = e;
        REQUIRE(pid.compute(vector_t::Constant(e), vector_t::Zero(), Scalar{1} / Scalar{100}).has_value());
        const bool oscillating = pid.oscillating();
        if(oscillating && run.first_detection < 0)
            run.first_detection = k;
        if(run.first_detection >= 0 && !oscillating)
            run.held_after_first_detection = false;
    }
    run.final_metric = static_cast<double>(pid.template metric<ctrlpp::oscillation_detect>()[0]);
    return run;
}

// Hagglund's x right after the n-th load event is the sum over the loads so far of gamma^(k_n - k_i), with
// gamma = 1 - h/T_sup = 1 - 2e-5 per step. Half-period spacing (500 steps) gives gamma^500 ~ 0.990, so ten
// loads sum to about 9.56 < 10 = n_lim and eleven to about 10.46 > 10: the verdict first holds at the
// eleventh load. The sum is formed here from the observed load steps rather than the nominal spacing.
auto load_rate_after(const std::vector<int>& loads, std::size_t n) -> double
{
    const double gamma = 1.0 - 0.01 / 500.0;
    double x = 0.0;
    for(std::size_t i = 0; i < n; ++i)
        x += std::pow(gamma, loads[n - 1] - loads[i]);
    return x;
}

void require_detection_at_eleventh_load(const detector_run& run)
{
    REQUIRE(run.sign_changes.size() >= 11);
    REQUIRE(load_rate_after(run.sign_changes, 10) < 10.0);
    REQUIRE(load_rate_after(run.sign_changes, 11) > 10.0);
    CAPTURE(run.sign_changes[9], run.sign_changes[10], run.first_detection);
    REQUIRE(run.first_detection == run.sign_changes[10]);
    REQUIRE(run.held_after_first_detection);
}

}

TEST_CASE("IAE accumulates integral of |error| * dt",
    "[pid][siso][perf-assessment][iae]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::IAE>;
    using PaPid = ctrlpp::pid<double, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    auto pid = ctrlpp::test::constructed(PaPid::create(cfg));

    // Constant error=1.0, 10 steps of dt=0.1
    // IAE = sum(|1.0| * 0.1) = 1.0
    for (int i = 0; i < 10; ++i)
        REQUIRE(pid.compute(vec1(1.0), vec1(0.0), 0.1).has_value());

    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[0], WithinAbs(1.0, 1e-10));
}

TEST_CASE("ISE accumulates integral of error^2 * dt",
    "[pid][siso][perf-assessment][ise]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::ISE>;
    using PaPid = ctrlpp::pid<double, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    auto pid = ctrlpp::test::constructed(PaPid::create(cfg));

    // Constant error=2.0, 10 steps of dt=0.1
    // ISE = sum(4.0 * 0.1) = 4.0
    for (int i = 0; i < 10; ++i)
        REQUIRE(pid.compute(vec1(2.0), vec1(0.0), 0.1).has_value());

    REQUIRE_THAT(pid.metric<ctrlpp::ISE>()[0], WithinAbs(4.0, 1e-10));
}

TEST_CASE("ITAE accumulates integral of t * |error| * dt",
    "[pid][siso][perf-assessment][itae]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::ITAE>;
    using PaPid = ctrlpp::pid<double, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    auto pid = ctrlpp::test::constructed(PaPid::create(cfg));

    // Constant error=1.0, dt=0.1
    // After step k: accumulated_time = k * 0.1
    // ITAE = sum(t_k * |1.0| * 0.1)
    // t_k = 0.1, 0.2, ..., 1.0 (accumulated AFTER dt update)
    // ITAE = 0.1*(0.1 + 0.2 + ... + 1.0) = 0.1 * 5.5 = 0.55
    for (int i = 0; i < 10; ++i)
        REQUIRE(pid.compute(vec1(1.0), vec1(0.0), 0.1).has_value());

    REQUIRE_THAT(pid.metric<ctrlpp::ITAE>()[0], WithinAbs(0.55, 1e-10));
}

TEST_CASE("Multiple metrics accumulate simultaneously",
    "[pid][siso][perf-assessment][multi]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::IAE, ctrlpp::ISE>;
    using PaPid = ctrlpp::pid<double, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    auto pid = ctrlpp::test::constructed(PaPid::create(cfg));

    for (int i = 0; i < 10; ++i)
        REQUIRE(pid.compute(vec1(2.0), vec1(0.0), 0.1).has_value());

    // IAE = sum(|2.0| * 0.1) = 2.0
    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[0], WithinAbs(2.0, 1e-10));
    // ISE = sum(4.0 * 0.1) = 4.0
    REQUIRE_THAT(pid.metric<ctrlpp::ISE>()[0], WithinAbs(4.0, 1e-10));
}

// A 10 s limit cycle of amplitude 1: each half period's IAE is about 2A/omega_u = 10/pi, twice
// IAE_lim = 2a/omega_u = 5/pi, so every sign change is a load event.
TEMPLATE_TEST_CASE("oscillation_detect flags a 10 s limit cycle at the eleventh load event",
    "[pid][siso][perf-assessment][oscillation]", float, double)
{
    using Scalar = TestType;
    auto const sine = [](int k) {
        return std::sin(Scalar{2} * std::numbers::pi_v<Scalar> * static_cast<Scalar>(k % steps_per_period) / Scalar{1000});
    };
    auto const run = run_detector<Scalar>(Scalar{1} / Scalar{2}, 10 * steps_per_period, sine);
    require_detection_at_eleventh_load(run);
}

// At amplitude a/2 the half-period IAE is IAE_lim/2. Sampling moves it by the (omega h)^2/12 ~ 3.3e-6
// relative discretization effect plus at most one sample, A h = 2.5e-3, at the crossings; both are far
// below the factor of two, so no sign change is a load event and x stays exactly zero.
TEMPLATE_TEST_CASE("oscillation_detect ignores a limit cycle at half the acceptable amplitude",
    "[pid][siso][perf-assessment][oscillation]", float, double)
{
    using Scalar = TestType;
    auto const sine = [](int k) {
        return std::sin(Scalar{2} * std::numbers::pi_v<Scalar> * static_cast<Scalar>(k % steps_per_period) / Scalar{1000})
            / Scalar{4};
    };
    auto const run = run_detector<Scalar>(Scalar{1} / Scalar{2}, 10 * steps_per_period, sine);
    REQUIRE(run.sign_changes.size() >= 19);
    REQUIRE(run.first_detection == -1);
    REQUIRE(run.final_metric == 0.0);
}

// 1e-6 noise of random sign at 100 Hz for 1000 s: the IAE of the whole run is 1e5 * 1e-6 * 1e-2 = 1e-3,
// below IAE_lim = 2a/omega_u ~ 3.2e-2 for a = 1/100, so no interval between sign changes can be a load.
TEMPLATE_TEST_CASE("oscillation_detect never flags 1e-6 noise at 100 Hz",
    "[pid][siso][perf-assessment][oscillation]", float, double)
{
    using Scalar = TestType;
    std::mt19937 bits{20261011u};
    auto const noise = [&bits](int) { return (bits() & 1u) != 0u ? static_cast<Scalar>(1e-6) : -static_cast<Scalar>(1e-6); };
    auto const run = run_detector<Scalar>(Scalar{1} / Scalar{100}, 100 * steps_per_period, noise);
    REQUIRE(run.sign_changes.size() > 1000);
    REQUIRE(run.first_detection == -1);
    REQUIRE(run.final_metric == 0.0);
}

TEMPLATE_TEST_CASE("oscillation_detect flags a slow square wave with the noise case's settings",
    "[pid][siso][perf-assessment][oscillation]", float, double)
{
    using Scalar = TestType;
    auto const square = [](int k) { return k % steps_per_period < steps_per_period / 2 ? Scalar{1} : Scalar{-1}; };
    auto const run = run_detector<Scalar>(Scalar{1} / Scalar{100}, 10 * steps_per_period, square);
    require_detection_at_eleventh_load(run);
}

TEST_CASE("oscillation_detect: constant error sign -> not oscillating",
    "[pid][siso][perf-assessment][oscillation]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::oscillation_detect>;
    using PaPid = ctrlpp::pid<double, 1, PA>;
    auto pid = ctrlpp::test::constructed(PaPid::create(detector_config<double, PA>(1.0 / 32.0, 1.0)));

    for (int i = 0; i < 10; ++i)
        REQUIRE(pid.compute(vec1(1.0), vec1(0.0), 0.1).has_value());

    REQUIRE(pid.metric<ctrlpp::oscillation_detect>()[0] == 0.0);
    REQUIRE(pid.oscillating() == false);
}

// Alternating +-1 at dt = 0.1: every interval between sign changes holds IAE 0.1 > IAE_lim = 2a/omega_u
// = 1/16, so each of the 9 sign changes is a load and x = sum of gamma^j, j = 0..8, lies in (9 gamma^8, 9].
TEST_CASE("oscillation_detect + IAE both accumulate",
    "[pid][siso][perf-assessment][oscillation][iae]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::oscillation_detect, ctrlpp::IAE>;
    using PaPid = ctrlpp::pid<double, 1, PA>;
    auto pid = ctrlpp::test::constructed(PaPid::create(detector_config<double, PA>(1.0 / 32.0, 1.0)));

    for (int i = 0; i < 10; ++i) {
        double sp = (i % 2 == 0) ? 1.0 : -1.0;
        REQUIRE(pid.compute(vec1(sp), vec1(0.0), 0.1).has_value());
    }

    const double gamma = 1.0 - 0.1 / (50.0 * 2.0 * std::numbers::pi);
    REQUIRE(pid.metric<ctrlpp::oscillation_detect>()[0] > 9.0 * std::pow(gamma, 8));
    REQUIRE(pid.metric<ctrlpp::oscillation_detect>()[0] <= 9.0);
    REQUIRE(pid.oscillating() == false);
    // IAE = sum(|1.0| * 0.1) = 1.0
    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[0], WithinAbs(1.0, 1e-10));
}

TEST_CASE("reset_metrics clears all metric accumulators",
    "[pid][siso][perf-assessment][reset]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::IAE, ctrlpp::oscillation_detect>;
    using PaPid = ctrlpp::pid<double, 1, PA>;
    auto pid = ctrlpp::test::constructed(PaPid::create(detector_config<double, PA>(1.0 / 32.0, 1.0)));

    for (int i = 0; i < 10; ++i) {
        double sp = (i % 2 == 0) ? 1.0 : -1.0;
        REQUIRE(pid.compute(vec1(sp), vec1(0.0), 0.1).has_value());
    }

    REQUIRE(pid.metric<ctrlpp::IAE>()[0] > 0.0);
    REQUIRE(pid.metric<ctrlpp::oscillation_detect>()[0] > 0.0);

    pid.reset_metrics();

    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[0], WithinAbs(0.0, tol));
    REQUIRE(pid.metric<ctrlpp::oscillation_detect>()[0] == 0.0);
    REQUIRE(pid.oscillating() == false);
}

TEST_CASE("IAE accumulates absolute value for negative errors",
    "[pid][siso][perf-assessment][iae][negative]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::IAE>;
    using PaPid = ctrlpp::pid<double, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    auto pid = ctrlpp::test::constructed(PaPid::create(cfg));

    // Negative error: sp=0, meas=1 -> e=-1
    for (int i = 0; i < 10; ++i)
        REQUIRE(pid.compute(vec1(0.0), vec1(1.0), 0.1).has_value());

    // IAE = sum(|-1| * 0.1) = 1.0 (same as positive error)
    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[0], WithinAbs(1.0, 1e-10));
}

TEST_CASE("ITAE weights later errors more heavily than earlier errors",
    "[pid][siso][perf-assessment][itae][negative]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::ITAE>;
    using PaPid = ctrlpp::pid<double, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    auto pid = ctrlpp::test::constructed(PaPid::create(cfg));

    // Negative error: sp=0, meas=2 -> e=-2
    // ITAE = sum(t_k * |e| * dt), |e|=2
    // t_k = 0.1, 0.2, ..., 1.0
    // ITAE = 2 * 0.1 * (0.1+0.2+...+1.0) = 2 * 0.1 * 5.5 = 1.1
    for (int i = 0; i < 10; ++i)
        REQUIRE(pid.compute(vec1(0.0), vec1(2.0), 0.1).has_value());

    REQUIRE_THAT(pid.metric<ctrlpp::ITAE>()[0], WithinAbs(1.1, 1e-10));
}

// +1, 0, -1, 0, +1 with every crossing a load (IAE 0.1 > 1/16): an exact zero neither crosses nor
// restarts the interval, so exactly two loads arrive and x = 1 + gamma^2 lies strictly between 1 and 2.
// Treating a zero as a sign change would add loads and push x past 2.
TEST_CASE("oscillation_detect skips zero-error steps (no false crossing)",
    "[pid][siso][perf-assessment][oscillation][zero-error]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::oscillation_detect>;
    using PaPid = ctrlpp::pid<double, 1, PA>;
    auto pid = ctrlpp::test::constructed(PaPid::create(detector_config<double, PA>(1.0 / 32.0, 1.0)));

    REQUIRE(pid.compute(vec1(1.0), vec1(0.0), 0.1).has_value());
    REQUIRE(pid.compute(vec1(0.0), vec1(0.0), 0.1).has_value());
    REQUIRE(pid.compute(vec1(0.0), vec1(1.0), 0.1).has_value());
    REQUIRE(pid.compute(vec1(0.0), vec1(0.0), 0.1).has_value());
    REQUIRE(pid.compute(vec1(1.0), vec1(0.0), 0.1).has_value());

    REQUIRE(pid.metric<ctrlpp::oscillation_detect>()[0] > 1.0);
    REQUIRE(pid.metric<ctrlpp::oscillation_detect>()[0] < 2.0);
}

TEST_CASE("oscillating() is false before any cycle",
    "[pid][siso][perf-assessment][oscillation][edge-case]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::oscillation_detect>;
    using PaPid = ctrlpp::pid<double, 1, PA>;
    auto pid = ctrlpp::test::constructed(PaPid::create(detector_config<double, PA>(1.0 / 32.0, 1.0)));

    REQUIRE(pid.metric<ctrlpp::oscillation_detect>()[0] == 0.0);
    REQUIRE(pid.oscillating() == false);
}

TEST_CASE("All four metrics accumulate simultaneously",
    "[pid][siso][perf-assessment][all-metrics]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::IAE, ctrlpp::ISE, ctrlpp::ITAE, ctrlpp::oscillation_detect>;
    using PaPid = ctrlpp::pid<double, 1, PA>;
    auto pid = ctrlpp::test::constructed(PaPid::create(detector_config<double, PA>(1.0 / 32.0, 1.0)));

    // Alternating error: 1, -1, 1 at dt=0.1
    REQUIRE(pid.compute(vec1(1.0), vec1(0.0), 0.1).has_value());  // e=1, t=0.1
    REQUIRE(pid.compute(vec1(0.0), vec1(1.0), 0.1).has_value());  // e=-1, t=0.2
    REQUIRE(pid.compute(vec1(1.0), vec1(0.0), 0.1).has_value());  // e=1, t=0.3

    // IAE = |1|*0.1 + |-1|*0.1 + |1|*0.1 = 0.3
    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[0], WithinAbs(0.3, 1e-10));
    // ISE = 1*0.1 + 1*0.1 + 1*0.1 = 0.3
    REQUIRE_THAT(pid.metric<ctrlpp::ISE>()[0], WithinAbs(0.3, 1e-10));
    // ITAE = 0.1*1*0.1 + 0.2*1*0.1 + 0.3*1*0.1 = 0.01+0.02+0.03 = 0.06
    REQUIRE_THAT(pid.metric<ctrlpp::ITAE>()[0], WithinAbs(0.06, 1e-10));
    // Two sign changes, each a load (IAE 0.1 > 1/16): x = gamma + 1, strictly between 1 and 2.
    REQUIRE(pid.metric<ctrlpp::oscillation_detect>()[0] > 1.0);
    REQUIRE(pid.metric<ctrlpp::oscillation_detect>()[0] < 2.0);
}

TEST_CASE("reset() clears performance metrics when perf_assessment is enabled",
    "[pid][siso][perf-assessment][reset-full]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::IAE, ctrlpp::ISE, ctrlpp::ITAE>;
    using PaPid = ctrlpp::pid<double, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    auto pid = ctrlpp::test::constructed(PaPid::create(cfg));

    for (int i = 0; i < 5; ++i)
        REQUIRE(pid.compute(vec1(1.0), vec1(0.0), 0.1).has_value());

    REQUIRE(pid.metric<ctrlpp::IAE>()[0] > 0.0);
    REQUIRE(pid.metric<ctrlpp::ISE>()[0] > 0.0);
    REQUIRE(pid.metric<ctrlpp::ITAE>()[0] > 0.0);

    pid.reset();  // full reset (not just reset_metrics)

    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[0], WithinAbs(0.0, tol));
    REQUIRE_THAT(pid.metric<ctrlpp::ISE>()[0], WithinAbs(0.0, tol));
    REQUIRE_THAT(pid.metric<ctrlpp::ITAE>()[0], WithinAbs(0.0, tol));
}
