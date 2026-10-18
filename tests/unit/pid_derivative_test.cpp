#include "hardening_helpers.h"
#include "ctrlpp/pid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <type_traits>

using Catch::Matchers::WithinAbs;

namespace {

using SisoPid = ctrlpp::pid<double, 1>;
using Vec1 = ctrlpp::Vector<double, 1>;

constexpr double Ts = 0.01;
constexpr double tol = 1e-12;

Vec1 vec1(double v) { Vec1 r; r << v; return r; }

template <typename Scalar>
using filtered_pid = ctrlpp::pid<Scalar, 1, ctrlpp::deriv_filter>;

template <typename Scalar>
auto filtered_config(Scalar kp, Scalar kd, Scalar n, std::optional<Scalar> tf) -> typename filtered_pid<Scalar>::config_type
{
    using DfPid = filtered_pid<Scalar>;
    typename DfPid::config_type cfg{};
    cfg.kp = DfPid::vector_t::Constant(kp);
    cfg.kd = DfPid::vector_t::Constant(kd);
    cfg.template policy<ctrlpp::deriv_filter>().n = {n};
    cfg.template policy<ctrlpp::deriv_filter>().tf = {tf};
    return cfg;
}

template <typename Scalar>
auto filter_refusal(Scalar kp, Scalar kd, Scalar n, std::optional<Scalar> tf) -> ctrlpp::pid_config_error
{
    const auto created = filtered_pid<Scalar>::create(filtered_config(kp, kd, n, tf));
    REQUIRE_FALSE(created.has_value());
    return created.error();
}

}

TEST_CASE("PD derivative on measurement (default)", "[pid][siso][derivative]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.kd = vec1(0.1);
    auto pid = ctrlpp::test::constructed(SisoPid::create(cfg));

    // Step 1: first step, derivative is zero
    auto u1 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    // P=1.0, D=0 (first step)
    REQUIRE_THAT(u1[0], WithinAbs(1.0, tol));

    // Step 2: measurement changes from 0.0 to 0.2
    // D = -Kd * (meas_new - meas_old) / dt = -0.1 * (0.2 - 0.0) / 0.01 = -2.0
    auto u2 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.2), Ts));
    // P = 1.0 * (1.0 - 0.2) = 0.8
    REQUIRE_THAT(u2[0], WithinAbs(0.8 - 2.0, tol));
}

TEST_CASE("PD derivative on error", "[pid][siso][derivative]")
{
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.kd = vec1(0.1);
    cfg.derivative_on_error = true;
    auto pid = ctrlpp::test::constructed(SisoPid::create(cfg));

    // Step 1: first step, derivative is zero
    auto u1 = ctrlpp::test::commanded(pid.compute(vec1(1.0), vec1(0.0), Ts));
    REQUIRE_THAT(u1[0], WithinAbs(1.0, tol));

    // Step 2: setpoint jumps from 1.0 to 2.0, meas stays at 0.0
    // e_new = 2.0 - 0.0 = 2.0, e_old = 1.0
    // D = Kd * (e_new - e_old) / dt = 0.1 * (2.0 - 1.0) / 0.01 = 10.0
    auto u2 = ctrlpp::test::commanded(pid.compute(vec1(2.0), vec1(0.0), Ts));
    // P = 1.0 * 2.0 = 2.0
    REQUIRE_THAT(u2[0], WithinAbs(2.0 + 10.0, tol));
}

TEST_CASE("deriv_filter reduces peak derivative on step input",
    "[pid][siso][deriv-filter]")
{
    // Unfiltered PID
    SisoPid::config_type cfg_unfiltered{};
    cfg_unfiltered.kp = vec1(1.0);
    cfg_unfiltered.kd = vec1(1.0);
    auto pid_unfiltered = ctrlpp::test::constructed(SisoPid::create(cfg_unfiltered));

    // Filtered PID
    using DfPid = ctrlpp::pid<double, 1, ctrlpp::deriv_filter>;
    DfPid::config_type cfg_filtered{};
    cfg_filtered.kp = vec1(1.0);
    cfg_filtered.kd = vec1(1.0);
    cfg_filtered.template policy<ctrlpp::deriv_filter>().n = {10.0};
    auto pid_filtered = ctrlpp::test::constructed(DfPid::create(cfg_filtered));

    // Step 1: first step, D=0 for both
    REQUIRE(pid_unfiltered.compute(vec1(0.0), vec1(0.0), Ts).has_value());
    REQUIRE(pid_filtered.compute(vec1(0.0), vec1(0.0), Ts).has_value());

    // Step 2: measurement step from 0 to 1 -> large derivative
    auto u_unfiltered = ctrlpp::test::commanded(pid_unfiltered.compute(vec1(0.0), vec1(1.0), Ts));
    auto u_filtered = ctrlpp::test::commanded(pid_filtered.compute(vec1(0.0), vec1(1.0), Ts));

    // Unfiltered: D = -Kd * (1-0)/dt = -1/0.01 = -100
    // P = 1.0 * (0-1) = -1
    // u_unfiltered = -1 + (-100) = -101
    REQUIRE_THAT(u_unfiltered[0], WithinAbs(-101.0, tol));

    // Filtered: D should have smaller magnitude (filter dampens the spike)
    // Tf = Kd / (Kp * N) = 1.0 / (1.0 * 10.0) = 0.1
    // alpha = Tf / (Tf + Ts) = 0.1 / 0.11 ~ 0.9091
    // D_filtered = alpha * 0 + (1 - alpha) * (-100) ~ -9.09
    double tf = 1.0 / (1.0 * 10.0);
    double alpha = tf / (tf + Ts);
    double d_filtered = (1.0 - alpha) * (-100.0);
    // P = -1
    REQUIRE_THAT(u_filtered[0], WithinAbs(-1.0 + d_filtered, 1e-6));

    // Peak derivative magnitude is smaller with filter
    REQUIRE(std::abs(u_filtered[0]) < std::abs(u_unfiltered[0]));
}

TEST_CASE("deriv_filter with N parameter smooths derivative over multiple steps",
    "[pid][siso][deriv-filter]")
{
    using DfPid = ctrlpp::pid<double, 1, ctrlpp::deriv_filter>;
    DfPid::config_type cfg{};
    cfg.kp = vec1(2.0);
    cfg.kd = vec1(0.5);
    cfg.template policy<ctrlpp::deriv_filter>().n = {20.0};
    auto pid = ctrlpp::test::constructed(DfPid::create(cfg));

    // Tf = Kd / (Kp * N) = 0.5 / (2.0 * 20.0) = 0.0125
    double tf_val = 0.5 / (2.0 * 20.0);
    double alpha = tf_val / (tf_val + Ts);

    // Step 1: D=0 (first step)
    REQUIRE(pid.compute(vec1(0.0), vec1(0.0), Ts).has_value());

    // Step 2: measurement step to 1.0
    // D_raw = -Kd * (1-0)/dt = -0.5 / 0.01 = -50
    // D_filt = alpha*0 + (1-alpha)*(-50)
    double d_filt = (1.0 - alpha) * (-50.0);
    auto u2 = ctrlpp::test::commanded(pid.compute(vec1(0.0), vec1(1.0), Ts));
    double p2 = 2.0 * (0.0 - 1.0);
    REQUIRE_THAT(u2[0], WithinAbs(p2 + d_filt, 1e-6));

    // Step 3: measurement stays at 1.0
    // D_raw = -Kd * (1-1)/dt = 0
    // D_filt = alpha*d_filt + (1-alpha)*0
    double d_filt2 = alpha * d_filt;
    auto u3 = ctrlpp::test::commanded(pid.compute(vec1(0.0), vec1(1.0), Ts));
    double p3 = 2.0 * (0.0 - 1.0);
    REQUIRE_THAT(u3[0], WithinAbs(p3 + d_filt2, 1e-6));
}

TEST_CASE("Without deriv_filter policy, unfiltered derivative is used",
    "[pid][siso][deriv-filter][compile]")
{
    // This ensures that without deriv_filter, no filter overhead occurs
    SisoPid::config_type cfg{};
    cfg.kp = vec1(1.0);
    cfg.kd = vec1(1.0);
    auto pid = ctrlpp::test::constructed(SisoPid::create(cfg));

    REQUIRE(pid.compute(vec1(0.0), vec1(0.0), Ts).has_value());
    auto u = ctrlpp::test::commanded(pid.compute(vec1(0.0), vec1(1.0), Ts));
    // D = -1.0 * (1-0)/0.01 = -100
    // P = 1.0 * (0-1) = -1
    REQUIRE_THAT(u[0], WithinAbs(-101.0, tol));
}

TEMPLATE_TEST_CASE("deriv_filter refuses a derivative channel whose filter time constant is undefined",
    "[pid][siso][deriv-filter]", float, double)
{
    using Scalar = TestType;
    using DfPid = filtered_pid<Scalar>;
    const Scalar kd = static_cast<Scalar>(0.1);
    const auto none = std::optional<Scalar>{};
    REQUIRE(filter_refusal<Scalar>(Scalar{0}, kd, Scalar{10}, none) == ctrlpp::pid_config_error::undefined_derivative_filter);
    REQUIRE(filter_refusal<Scalar>(Scalar{1}, kd, Scalar{0}, none) == ctrlpp::pid_config_error::undefined_derivative_filter);
    REQUIRE(filter_refusal<Scalar>(Scalar{1}, -kd, Scalar{10}, none) == ctrlpp::pid_config_error::undefined_derivative_filter);

    for(const Scalar n : {Scalar{-1}, Scalar{0}, Scalar{10}})
        REQUIRE(DfPid::create(filtered_config<Scalar>(Scalar{1}, Scalar{0}, n, none)).has_value());

    for(const Scalar tf : {Scalar{0}, Scalar{-1}, std::numeric_limits<Scalar>::quiet_NaN(), std::numeric_limits<Scalar>::infinity()})
        REQUIRE(filter_refusal<Scalar>(Scalar{1}, kd, Scalar{10}, tf) == ctrlpp::pid_config_error::invalid_derivative_filter);
}

TEMPLATE_TEST_CASE("deriv_filter with an explicit time constant filters a channel without proportional gain",
    "[pid][siso][deriv-filter]", float, double)
{
    using Scalar = TestType;
    using DfPid = filtered_pid<Scalar>;
    using vector_t = typename DfPid::vector_t;
    const Scalar kd = static_cast<Scalar>(0.1);
    const Scalar tf = static_cast<Scalar>(0.02);
    const Scalar dt = static_cast<Scalar>(0.01);
    auto controller = ctrlpp::test::constructed(DfPid::create(filtered_config<Scalar>(Scalar{0}, kd, Scalar{0}, tf)));
    REQUIRE(controller.compute(vector_t::Zero(), vector_t::Zero(), dt).has_value());
    const auto u = ctrlpp::test::commanded(controller.compute(vector_t::Zero(), vector_t::Constant(Scalar{1}), dt));

    const Scalar alpha = tf / (tf + dt);
    const Scalar raw = -(kd * ((Scalar{1} - Scalar{0}) / dt));
    REQUIRE(u[0] == alpha * Scalar{0} + (Scalar{1} - alpha) * raw);
}

// The default time constant is Td / N = kd / (kp * n). Here kp * n = 20 is exact and the quotient rounds once, so a
// controller handed that quotient runs the defaulted one's arithmetic on any platform and must agree with it bit for bit.
TEMPLATE_TEST_CASE("deriv_filter omitting the time constant is bit-identical to naming kd / (kp * n)", "[pid][siso][deriv-filter]", float, double)
{
    using Scalar = TestType;
    using DfPid = filtered_pid<Scalar>;
    using vector_t = typename DfPid::vector_t;
    const Scalar kp = 2;
    const Scalar kd = static_cast<Scalar>(0.1);
    const Scalar n = 10;
    auto defaulted = ctrlpp::test::constructed(DfPid::create(filtered_config<Scalar>(kp, kd, n, std::nullopt)));
    auto named = ctrlpp::test::constructed(DfPid::create(filtered_config<Scalar>(kp, kd, n, kd / (kp * n))));
    constexpr std::array<double, 8> measurements{0, 1, 1, 0.5, 0.25, 0.75, 0.75, -0.5};
    for(std::size_t k = 0; k < measurements.size(); ++k)
    {
        CAPTURE(k);
        const auto meas = vector_t::Constant(static_cast<Scalar>(measurements[k]));
        const auto u_defaulted = ctrlpp::test::commanded(defaulted.compute(vector_t::Zero(), meas, static_cast<Scalar>(0.01)));
        const auto u_named = ctrlpp::test::commanded(named.compute(vector_t::Zero(), meas, static_cast<Scalar>(0.01)));
        REQUIRE(u_defaulted[0] == u_named[0]);
    }
}
