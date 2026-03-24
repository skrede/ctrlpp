#include <ctrlpp/pid.h>
#include <ctrlpp/pid_policies.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;

namespace
{

constexpr double dt = 0.01;
constexpr double tol = 1e-10;

} // namespace

TEST_CASE("SISO PID with Eigen - P-only", "[pid][eigen][siso]")
{
    using pid = ctrlpp::pid<double, 1, 1, 1>;
    using Vec = pid::vector_t;

    pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.5);
    pid ctrl(cfg);

    auto u = ctrl.compute(Vec::Constant(1.0), Vec::Constant(0.0), dt);
    REQUIRE_THAT(u[0], WithinAbs(2.5, tol));
}

TEST_CASE("SISO PID with Eigen - PI", "[pid][eigen][siso]")
{
    using pid = ctrlpp::pid<double, 1, 1, 1>;
    using Vec = pid::vector_t;

    pid::config_type cfg{};
    cfg.kp = Vec::Constant(1.0);
    cfg.ki = Vec::Constant(0.5);
    pid ctrl(cfg);

    auto u1 = ctrl.compute(Vec::Constant(1.0), Vec::Constant(0.0), dt);
    REQUIRE_THAT(u1[0], WithinAbs(1.005, tol));

    auto u2 = ctrl.compute(Vec::Constant(1.0), Vec::Constant(0.0), dt);
    REQUIRE_THAT(u2[0], WithinAbs(1.01, tol));
}

TEST_CASE("SISO PID with Eigen - PID", "[pid][eigen][siso]")
{
    using pid = ctrlpp::pid<double, 1, 1, 1>;
    using Vec = pid::vector_t;

    pid::config_type cfg{};
    cfg.kp = Vec::Constant(1.0);
    cfg.ki = Vec::Constant(0.5);
    cfg.kd = Vec::Constant(0.1);
    cfg.derivative_on_error = true;
    pid ctrl(cfg);

    auto u1 = ctrl.compute(Vec::Constant(1.0), Vec::Constant(0.0), dt);
    REQUIRE_THAT(u1[0], WithinAbs(1.005, tol));

    auto u2 = ctrl.compute(Vec::Constant(1.0), Vec::Constant(0.0), dt);
    REQUIRE_THAT(u2[0], WithinAbs(1.01, tol));
}

TEST_CASE("MIMO PID with Eigen - 2-channel independence", "[pid][eigen][mimo]")
{
    using pid = ctrlpp::pid<double, 2, 2, 2>;
    using Vec = pid::vector_t;

    pid::config_type cfg{};
    cfg.kp = (Vec() << 1.0, 2.0).finished();
    cfg.ki = (Vec() << 0.0, 0.0).finished();
    cfg.kd = (Vec() << 0.0, 0.0).finished();

    pid ctrl(cfg);

    Vec sp = (Vec() << 1.0, 3.0).finished();
    Vec meas = (Vec() << 0.0, 1.0).finished();

    auto u = ctrl.compute(sp, meas, dt);

    REQUIRE_THAT(u[0], WithinAbs(1.0, tol));
    REQUIRE_THAT(u[1], WithinAbs(4.0, tol));
}

TEST_CASE("MIMO PID with Eigen - per-channel PI gains", "[pid][eigen][mimo]")
{
    using pid = ctrlpp::pid<double, 2, 2, 2>;
    using Vec = pid::vector_t;

    pid::config_type cfg{};
    cfg.kp = (Vec() << 1.0, 2.0).finished();
    cfg.ki = (Vec() << 0.5, 1.0).finished();

    pid ctrl(cfg);

    Vec sp = (Vec() << 1.0, 1.0).finished();
    Vec meas = Vec::Zero();

    auto u = ctrl.compute(sp, meas, dt);

    REQUIRE_THAT(u[0], WithinAbs(1.005, tol));
    REQUIRE_THAT(u[1], WithinAbs(2.01, tol));

    auto u2 = ctrl.compute(sp, meas, dt);
    REQUIRE_THAT(u2[0], WithinAbs(1.01, tol));
    REQUIRE_THAT(u2[1], WithinAbs(2.02, tol));
}

TEST_CASE("Full-featured PID with Eigen - anti_windup + deriv_filter + perf_assessment", "[pid][eigen][full]")
{
    using pid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::anti_windup<ctrlpp::back_calc>, ctrlpp::deriv_filter, ctrlpp::perf_assessment<ctrlpp::IAE>>;
    using Vec = pid::vector_t;

    pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.0);
    cfg.ki = Vec::Constant(1.0);
    cfg.kd = Vec::Constant(0.5);
    cfg.output_min = Vec::Constant(-10.0);
    cfg.output_max = Vec::Constant(10.0);
    cfg.template policy<ctrlpp::deriv_filter>().n = {10.0};

    pid ctrl(cfg);

    Vec sp = Vec::Constant(1.0);
    Vec meas = Vec::Constant(0.0);

    auto u1 = ctrl.compute(sp, meas, dt);
    REQUIRE(std::isfinite(u1[0]));
    REQUIRE(u1[0] > 0.0);

    auto u2 = ctrl.compute(sp, meas, dt);
    REQUIRE(std::isfinite(u2[0]));

    auto iae = ctrl.metric<ctrlpp::IAE>();
    REQUIRE(iae[0] > 0.0);
}

TEST_CASE("SISO with all composable policies - compile and run", "[pid][eigen][compose]")
{
    using pid =
        ctrlpp::pid<double, 1, 1, 1, ctrlpp::anti_windup<ctrlpp::back_calc>, ctrlpp::deriv_filter, ctrlpp::setpoint_filter, ctrlpp::pv_filter, ctrlpp::rate_limit, ctrlpp::perf_assessment<ctrlpp::IAE, ctrlpp::ISE>>;
    using Vec = pid::vector_t;

    pid::config_type cfg{};
    cfg.kp = Vec::Constant(1.0);
    cfg.ki = Vec::Constant(0.5);
    cfg.kd = Vec::Constant(0.1);
    cfg.output_min = Vec::Constant(-5.0);
    cfg.output_max = Vec::Constant(5.0);
    cfg.template policy<ctrlpp::deriv_filter>().n = {10.0};
    cfg.template policy<ctrlpp::setpoint_filter>().tf = {0.05};
    cfg.template policy<ctrlpp::pv_filter>().tf = {0.02};
    cfg.template policy<ctrlpp::rate_limit>().rate_max = {50.0};

    pid ctrl(cfg);

    Vec sp = Vec::Constant(1.0);
    Vec meas = Vec::Constant(0.0);

    for(int k = 0; k < 10; ++k)
    {
        auto u = ctrl.compute(sp, meas, dt);
        REQUIRE(std::isfinite(u[0]));
    }

    auto iae = ctrl.metric<ctrlpp::IAE>();
    auto ise = ctrl.metric<ctrlpp::ISE>();
    REQUIRE(iae[0] > 0.0);
    REQUIRE(ise[0] > 0.0);
}
