#include <ctrlpp/pid.h>
#include <ctrlpp/pid_policies.h>
#include <ctrlpp/eigen_linalg.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;

namespace {

using Policy = ctrlpp::EigenLinalgPolicy;
constexpr double dt = 0.01;
constexpr double tol = 1e-10;

}

TEST_CASE("SISO PID with Eigen - P-only", "[pid][eigen][siso]")
{
    using Pid = ctrlpp::Pid<Policy, double, 1, 1, 1>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.5);
    Pid pid(cfg);

    auto u = pid.compute(Vec::Constant(1.0), Vec::Constant(0.0), dt);
    REQUIRE_THAT(u[0], WithinAbs(2.5, tol));
}

TEST_CASE("SISO PID with Eigen - PI", "[pid][eigen][siso]")
{
    using Pid = ctrlpp::Pid<Policy, double, 1, 1, 1>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(1.0);
    cfg.ki = Vec::Constant(0.5);
    Pid pid(cfg);

    // Step 1: e=1.0, P=1.0, I=0.5*1.0*0.01=0.005
    auto u1 = pid.compute(Vec::Constant(1.0), Vec::Constant(0.0), dt);
    REQUIRE_THAT(u1[0], WithinAbs(1.005, tol));

    // Step 2: same error, integral accumulates
    auto u2 = pid.compute(Vec::Constant(1.0), Vec::Constant(0.0), dt);
    REQUIRE_THAT(u2[0], WithinAbs(1.01, tol));
}

TEST_CASE("SISO PID with Eigen - PID", "[pid][eigen][siso]")
{
    using Pid = ctrlpp::Pid<Policy, double, 1, 1, 1>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(1.0);
    cfg.ki = Vec::Constant(0.5);
    cfg.kd = Vec::Constant(0.1);
    cfg.derivative_on_error = true;
    Pid pid(cfg);

    // First step: D term is zero (first_step_ = true)
    auto u1 = pid.compute(Vec::Constant(1.0), Vec::Constant(0.0), dt);
    REQUIRE_THAT(u1[0], WithinAbs(1.005, tol));

    // Second step: e=1.0 still, de/dt=0 -> D=0
    auto u2 = pid.compute(Vec::Constant(1.0), Vec::Constant(0.0), dt);
    REQUIRE_THAT(u2[0], WithinAbs(1.01, tol));
}

TEST_CASE("MIMO PID with Eigen - 2-channel independence", "[pid][eigen][mimo]")
{
    using Pid = ctrlpp::Pid<Policy, double, 2, 2, 2>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = (Vec() << 1.0, 2.0).finished();
    cfg.ki = (Vec() << 0.0, 0.0).finished();
    cfg.kd = (Vec() << 0.0, 0.0).finished();

    Pid pid(cfg);

    Vec sp = (Vec() << 1.0, 3.0).finished();
    Vec meas = (Vec() << 0.0, 1.0).finished();

    auto u = pid.compute(sp, meas, dt);

    // Channel 0: kp=1.0 * e=1.0 = 1.0
    REQUIRE_THAT(u[0], WithinAbs(1.0, tol));
    // Channel 1: kp=2.0 * e=2.0 = 4.0
    REQUIRE_THAT(u[1], WithinAbs(4.0, tol));
}

TEST_CASE("MIMO PID with Eigen - per-channel PI gains", "[pid][eigen][mimo]")
{
    using Pid = ctrlpp::Pid<Policy, double, 2, 2, 2>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = (Vec() << 1.0, 2.0).finished();
    cfg.ki = (Vec() << 0.5, 1.0).finished();

    Pid pid(cfg);

    Vec sp = (Vec() << 1.0, 1.0).finished();
    Vec meas = Vec::Zero();

    auto u = pid.compute(sp, meas, dt);

    // Channel 0: P=1.0*1.0 + I=0.5*1.0*0.01 = 1.005
    REQUIRE_THAT(u[0], WithinAbs(1.005, tol));
    // Channel 1: P=2.0*1.0 + I=1.0*1.0*0.01 = 2.01
    REQUIRE_THAT(u[1], WithinAbs(2.01, tol));

    // Second step: integral accumulates independently
    auto u2 = pid.compute(sp, meas, dt);
    REQUIRE_THAT(u2[0], WithinAbs(1.01, tol));
    REQUIRE_THAT(u2[1], WithinAbs(2.02, tol));
}

TEST_CASE("Full-featured PID with Eigen - AntiWindup + DerivFilter + PerfAssessment",
          "[pid][eigen][full]")
{
    using Pid = ctrlpp::Pid<Policy, double, 1, 1, 1,
        ctrlpp::AntiWindup<ctrlpp::BackCalc>,
        ctrlpp::DerivFilter,
        ctrlpp::PerfAssessment<ctrlpp::IAE>>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.0);
    cfg.ki = Vec::Constant(1.0);
    cfg.kd = Vec::Constant(0.5);
    cfg.output_min = Vec::Constant(-10.0);
    cfg.output_max = Vec::Constant(10.0);
    cfg.template policy<ctrlpp::DerivFilter>().n = {10.0};

    Pid pid(cfg);

    // Run a few steps and verify output is finite and IAE accumulates
    Vec sp = Vec::Constant(1.0);
    Vec meas = Vec::Constant(0.0);

    auto u1 = pid.compute(sp, meas, dt);
    REQUIRE(std::isfinite(u1[0]));
    REQUIRE(u1[0] > 0.0);

    auto u2 = pid.compute(sp, meas, dt);
    REQUIRE(std::isfinite(u2[0]));

    auto iae = pid.metric<ctrlpp::IAE>();
    REQUIRE(iae[0] > 0.0);
}

TEST_CASE("SISO with all composable policies - compile and run", "[pid][eigen][compose]")
{
    using Pid = ctrlpp::Pid<Policy, double, 1, 1, 1,
        ctrlpp::AntiWindup<ctrlpp::BackCalc>,
        ctrlpp::DerivFilter,
        ctrlpp::SetpointFilter,
        ctrlpp::PvFilter,
        ctrlpp::RateLimit,
        ctrlpp::PerfAssessment<ctrlpp::IAE, ctrlpp::ISE>>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(1.0);
    cfg.ki = Vec::Constant(0.5);
    cfg.kd = Vec::Constant(0.1);
    cfg.output_min = Vec::Constant(-5.0);
    cfg.output_max = Vec::Constant(5.0);
    cfg.template policy<ctrlpp::DerivFilter>().n = {10.0};
    cfg.template policy<ctrlpp::SetpointFilter>().tf = {0.05};
    cfg.template policy<ctrlpp::PvFilter>().tf = {0.02};
    cfg.template policy<ctrlpp::RateLimit>().rate_max = {50.0};

    Pid pid(cfg);

    Vec sp = Vec::Constant(1.0);
    Vec meas = Vec::Constant(0.0);

    // Run 10 steps, verify finite output
    for (int k = 0; k < 10; ++k) {
        auto u = pid.compute(sp, meas, dt);
        REQUIRE(std::isfinite(u[0]));
    }

    auto iae = pid.metric<ctrlpp::IAE>();
    auto ise = pid.metric<ctrlpp::ISE>();
    REQUIRE(iae[0] > 0.0);
    REQUIRE(ise[0] > 0.0);
}
