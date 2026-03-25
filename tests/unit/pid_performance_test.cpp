#include "ctrlpp/pid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;

namespace {

using Vec1 = ctrlpp::Vector<double, 1>;

constexpr double tol = 1e-12;

Vec1 vec1(double v) { Vec1 r; r << v; return r; }

}

TEST_CASE("IAE accumulates integral of |error| * dt",
    "[pid][siso][perf-assessment][iae]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::IAE>;
    using PaPid = ctrlpp::pid<double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    PaPid pid(cfg);

    // Constant error=1.0, 10 steps of dt=0.1
    // IAE = sum(|1.0| * 0.1) = 1.0
    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(1.0), vec1(0.0), 0.1);

    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[0], WithinAbs(1.0, 1e-10));
}

TEST_CASE("ISE accumulates integral of error^2 * dt",
    "[pid][siso][perf-assessment][ise]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::ISE>;
    using PaPid = ctrlpp::pid<double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    PaPid pid(cfg);

    // Constant error=2.0, 10 steps of dt=0.1
    // ISE = sum(4.0 * 0.1) = 4.0
    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(2.0), vec1(0.0), 0.1);

    REQUIRE_THAT(pid.metric<ctrlpp::ISE>()[0], WithinAbs(4.0, 1e-10));
}

TEST_CASE("ITAE accumulates integral of t * |error| * dt",
    "[pid][siso][perf-assessment][itae]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::ITAE>;
    using PaPid = ctrlpp::pid<double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    PaPid pid(cfg);

    // Constant error=1.0, dt=0.1
    // After step k: accumulated_time = k * 0.1
    // ITAE = sum(t_k * |1.0| * 0.1)
    // t_k = 0.1, 0.2, ..., 1.0 (accumulated AFTER dt update)
    // ITAE = 0.1*(0.1 + 0.2 + ... + 1.0) = 0.1 * 5.5 = 0.55
    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(1.0), vec1(0.0), 0.1);

    REQUIRE_THAT(pid.metric<ctrlpp::ITAE>()[0], WithinAbs(0.55, 1e-10));
}

TEST_CASE("Multiple metrics accumulate simultaneously",
    "[pid][siso][perf-assessment][multi]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::IAE, ctrlpp::ISE>;
    using PaPid = ctrlpp::pid<double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    PaPid pid(cfg);

    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(2.0), vec1(0.0), 0.1);

    // IAE = sum(|2.0| * 0.1) = 2.0
    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[0], WithinAbs(2.0, 1e-10));
    // ISE = sum(4.0 * 0.1) = 4.0
    REQUIRE_THAT(pid.metric<ctrlpp::ISE>()[0], WithinAbs(4.0, 1e-10));
}

TEST_CASE("oscillation_detect counts zero-crossings and detects oscillation",
    "[pid][siso][perf-assessment][oscillation]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::oscillation_detect>;
    using PaPid = ctrlpp::pid<double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    PaPid pid(cfg);

    // Alternating error sign: +1, -1, +1, -1, ...  for 10 steps at dt=0.1
    // 9 zero-crossings over 1.0 second -> rate = 9/1.0 = 9 > threshold 5
    for (int i = 0; i < 10; ++i) {
        double sp = (i % 2 == 0) ? 1.0 : -1.0;
        pid.compute(vec1(sp), vec1(0.0), 0.1);
    }

    REQUIRE(pid.metric<ctrlpp::oscillation_detect>()[0] >= 9.0 - tol);
    REQUIRE(pid.oscillating() == true);
}

TEST_CASE("oscillation_detect: constant error sign -> not oscillating",
    "[pid][siso][perf-assessment][oscillation]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::oscillation_detect>;
    using PaPid = ctrlpp::pid<double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    PaPid pid(cfg);

    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(1.0), vec1(0.0), 0.1);

    REQUIRE_THAT(pid.metric<ctrlpp::oscillation_detect>()[0], WithinAbs(0.0, tol));
    REQUIRE(pid.oscillating() == false);
}

TEST_CASE("oscillation_detect + IAE both accumulate",
    "[pid][siso][perf-assessment][oscillation][iae]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::oscillation_detect, ctrlpp::IAE>;
    using PaPid = ctrlpp::pid<double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    PaPid pid(cfg);

    for (int i = 0; i < 10; ++i) {
        double sp = (i % 2 == 0) ? 1.0 : -1.0;
        pid.compute(vec1(sp), vec1(0.0), 0.1);
    }

    REQUIRE(pid.metric<ctrlpp::oscillation_detect>()[0] >= 9.0 - tol);
    // IAE = sum(|1.0| * 0.1) = 1.0
    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[0], WithinAbs(1.0, 1e-10));
}

TEST_CASE("reset_metrics clears all metric accumulators",
    "[pid][siso][perf-assessment][reset]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::IAE, ctrlpp::oscillation_detect>;
    using PaPid = ctrlpp::pid<double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    PaPid pid(cfg);

    for (int i = 0; i < 10; ++i) {
        double sp = (i % 2 == 0) ? 1.0 : -1.0;
        pid.compute(vec1(sp), vec1(0.0), 0.1);
    }

    REQUIRE(pid.metric<ctrlpp::IAE>()[0] > 0.0);
    REQUIRE(pid.metric<ctrlpp::oscillation_detect>()[0] > 0.0);

    pid.reset_metrics();

    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[0], WithinAbs(0.0, tol));
    REQUIRE_THAT(pid.metric<ctrlpp::oscillation_detect>()[0], WithinAbs(0.0, tol));
    REQUIRE(pid.oscillating() == false);
}

TEST_CASE("IAE accumulates absolute value for negative errors",
    "[pid][siso][perf-assessment][iae][negative]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::IAE>;
    using PaPid = ctrlpp::pid<double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    PaPid pid(cfg);

    // Negative error: sp=0, meas=1 -> e=-1
    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(0.0), vec1(1.0), 0.1);

    // IAE = sum(|-1| * 0.1) = 1.0 (same as positive error)
    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[0], WithinAbs(1.0, 1e-10));
}

TEST_CASE("ITAE weights later errors more heavily than earlier errors",
    "[pid][siso][perf-assessment][itae][negative]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::ITAE>;
    using PaPid = ctrlpp::pid<double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    PaPid pid(cfg);

    // Negative error: sp=0, meas=2 -> e=-2
    // ITAE = sum(t_k * |e| * dt), |e|=2
    // t_k = 0.1, 0.2, ..., 1.0
    // ITAE = 2 * 0.1 * (0.1+0.2+...+1.0) = 2 * 0.1 * 5.5 = 1.1
    for (int i = 0; i < 10; ++i)
        pid.compute(vec1(0.0), vec1(2.0), 0.1);

    REQUIRE_THAT(pid.metric<ctrlpp::ITAE>()[0], WithinAbs(1.1, 1e-10));
}

TEST_CASE("oscillation_detect skips zero-error steps (no false crossing)",
    "[pid][siso][perf-assessment][oscillation][zero-error]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::oscillation_detect>;
    using PaPid = ctrlpp::pid<double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    PaPid pid(cfg);

    // Sequence: +1, 0, -1, 0, +1
    // Zero errors should be ignored for crossing detection
    // Only +1 -> -1 and -1 -> +1 count (via non-zero signs)
    pid.compute(vec1(1.0), vec1(0.0), 0.1);   // e=+1, sign=+1
    pid.compute(vec1(0.0), vec1(0.0), 0.1);   // e=0, sign_e=0 -> prev_sign stays +1
    pid.compute(vec1(0.0), vec1(1.0), 0.1);   // e=-1, prev_sign=+1, crossing!
    pid.compute(vec1(0.0), vec1(0.0), 0.1);   // e=0, sign_e=0 -> prev_sign stays -1
    pid.compute(vec1(1.0), vec1(0.0), 0.1);   // e=+1, prev_sign=-1, crossing!

    REQUIRE_THAT(pid.metric<ctrlpp::oscillation_detect>()[0], WithinAbs(2.0, tol));
}

TEST_CASE("oscillating() returns false when accumulated_time is zero",
    "[pid][siso][perf-assessment][oscillation][edge-case]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::oscillation_detect>;
    using PaPid = ctrlpp::pid<double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    PaPid pid(cfg);

    // No steps computed yet -> accumulated_time = 0 -> oscillating() returns false
    REQUIRE(pid.oscillating() == false);
}

TEST_CASE("All four metrics accumulate simultaneously",
    "[pid][siso][perf-assessment][all-metrics]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::IAE, ctrlpp::ISE, ctrlpp::ITAE, ctrlpp::oscillation_detect>;
    using PaPid = ctrlpp::pid<double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    PaPid pid(cfg);

    // Alternating error: 1, -1, 1 at dt=0.1
    pid.compute(vec1(1.0), vec1(0.0), 0.1);  // e=1, t=0.1
    pid.compute(vec1(0.0), vec1(1.0), 0.1);  // e=-1, t=0.2
    pid.compute(vec1(1.0), vec1(0.0), 0.1);  // e=1, t=0.3

    // IAE = |1|*0.1 + |-1|*0.1 + |1|*0.1 = 0.3
    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[0], WithinAbs(0.3, 1e-10));
    // ISE = 1*0.1 + 1*0.1 + 1*0.1 = 0.3
    REQUIRE_THAT(pid.metric<ctrlpp::ISE>()[0], WithinAbs(0.3, 1e-10));
    // ITAE = 0.1*1*0.1 + 0.2*1*0.1 + 0.3*1*0.1 = 0.01+0.02+0.03 = 0.06
    REQUIRE_THAT(pid.metric<ctrlpp::ITAE>()[0], WithinAbs(0.06, 1e-10));
    // Oscillation crossings: +1 -> -1 (crossing), -1 -> +1 (crossing) = 2
    REQUIRE_THAT(pid.metric<ctrlpp::oscillation_detect>()[0], WithinAbs(2.0, tol));
}

TEST_CASE("reset() clears performance metrics when perf_assessment is enabled",
    "[pid][siso][perf-assessment][reset-full]")
{
    using PA = ctrlpp::perf_assessment<ctrlpp::IAE, ctrlpp::ISE, ctrlpp::ITAE>;
    using PaPid = ctrlpp::pid<double, 1, 1, 1, PA>;
    PaPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    PaPid pid(cfg);

    for (int i = 0; i < 5; ++i)
        pid.compute(vec1(1.0), vec1(0.0), 0.1);

    REQUIRE(pid.metric<ctrlpp::IAE>()[0] > 0.0);
    REQUIRE(pid.metric<ctrlpp::ISE>()[0] > 0.0);
    REQUIRE(pid.metric<ctrlpp::ITAE>()[0] > 0.0);

    pid.reset();  // full reset (not just reset_metrics)

    REQUIRE_THAT(pid.metric<ctrlpp::IAE>()[0], WithinAbs(0.0, tol));
    REQUIRE_THAT(pid.metric<ctrlpp::ISE>()[0], WithinAbs(0.0, tol));
    REQUIRE_THAT(pid.metric<ctrlpp::ITAE>()[0], WithinAbs(0.0, tol));
}
