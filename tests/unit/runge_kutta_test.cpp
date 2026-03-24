#include "ctrlpp/detail/runge_kutta.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

namespace
{

using Vec1 = ctrlpp::Vector<double, 1>;

// dx/dt = -x (exponential decay). Exact solution: x(t) = x0 * exp(-t)
struct exponential_decay
{
    auto operator()(const Vec1& x, const Vec1& /*u*/) const -> Vec1 { return -x; }
};

} // namespace

TEST_CASE("euler step on exponential decay")
{
    constexpr double dt = 0.01;
    Vec1 x;
    x << 1.0;
    Vec1 u;
    u << 0.0;

    auto x_next = ctrlpp::detail::euler_step(exponential_decay{}, x, u, dt);

    // Euler: x_next = x + dt * f(x, u) = 1 + 0.01 * (-1) = 0.99
    CHECK_THAT(x_next(0), Catch::Matchers::WithinAbs(0.99, 1e-14));
}

TEST_CASE("rk2 midpoint on exponential decay is more accurate than euler")
{
    constexpr double dt = 0.01;
    Vec1 x;
    x << 1.0;
    Vec1 u;
    u << 0.0;

    double exact = std::exp(-dt);

    auto x_euler = ctrlpp::detail::euler_step(exponential_decay{}, x, u, dt);
    auto x_rk2 = ctrlpp::detail::rk2_step(exponential_decay{}, x, u, dt);

    double err_euler = std::abs(x_euler(0) - exact);
    double err_rk2 = std::abs(x_rk2(0) - exact);

    CHECK(err_rk2 < err_euler);
}

TEST_CASE("rk4 on exponential decay matches exp(-dt) to high precision")
{
    constexpr double dt = 0.01;
    Vec1 x;
    x << 1.0;
    Vec1 u;
    u << 0.0;

    double exact = std::exp(-dt);
    auto x_rk4 = ctrlpp::detail::rk4_step(exponential_decay{}, x, u, dt);

    CHECK_THAT(x_rk4(0), Catch::Matchers::WithinAbs(exact, 1e-10));
}

TEST_CASE("order verification: halving dt reduces error by 2^order")
{
    Vec1 x;
    x << 1.0;
    Vec1 u;
    u << 0.0;

    auto compute_error = [&](auto step_fn, double dt)
    {
        auto x_next = step_fn(exponential_decay{}, x, u, dt);
        return std::abs(x_next(0) - std::exp(-dt));
    };

    // Use smaller dt so higher-order error terms are negligible
    constexpr double dt = 0.02;

    SECTION("euler is 1st order")
    {
        auto step = [](const auto& f, const Vec1& xx, const Vec1& uu, double h) { return ctrlpp::detail::euler_step(f, xx, uu, h); };
        double err1 = compute_error(step, dt);
        double err2 = compute_error(step, dt / 2.0);
        double ratio = err1 / err2;
        // 2^1 = 2; for smooth functions the effective ratio can exceed the
        // minimum due to error structure, so we verify the lower bound
        CHECK(ratio >= 2.0 - 0.3);
    }

    SECTION("rk2 is 2nd order")
    {
        auto step = [](const auto& f, const Vec1& xx, const Vec1& uu, double h) { return ctrlpp::detail::rk2_step(f, xx, uu, h); };
        double err1 = compute_error(step, dt);
        double err2 = compute_error(step, dt / 2.0);
        double ratio = err1 / err2;
        // 2^2 = 4 minimum
        CHECK(ratio >= 4.0 - 0.5);
    }

    SECTION("rk4 is 4th order")
    {
        auto step = [](const auto& f, const Vec1& xx, const Vec1& uu, double h) { return ctrlpp::detail::rk4_step(f, xx, uu, h); };
        double err1 = compute_error(step, dt);
        double err2 = compute_error(step, dt / 2.0);
        double ratio = err1 / err2;
        // 2^4 = 16 minimum
        CHECK(ratio >= 16.0 - 3.0);
    }
}

TEST_CASE("rk4 with input vector passthrough")
{
    // dx/dt = -x + u
    auto f = [](const Vec1& x, const Vec1& u) -> Vec1 { return Vec1{(-x + u).eval()}; };

    Vec1 x;
    x << 0.0;
    Vec1 u;
    u << 1.0;
    constexpr double dt = 0.01;

    auto x_next = ctrlpp::detail::rk4_step(f, x, u, dt);

    // dx/dt = -0 + 1 = 1 initially; x(dt) ~ dt for small dt
    CHECK(x_next(0) > 0.0);
    CHECK_THAT(x_next(0), Catch::Matchers::WithinAbs(dt, 1e-4));
}

TEST_CASE("rk45 adaptive integrates exponential decay over long interval")
{
    Vec1 x;
    x << 1.0;
    Vec1 u;
    u << 0.0;

    ctrlpp::detail::rk45_config<double> cfg{};
    cfg.atol = 1e-8;
    cfg.rtol = 1e-8;

    double t = 0.0;
    double dt = 0.1;
    constexpr double t_final = 1.0;

    while(t < t_final)
    {
        double dt_step = std::min(dt, t_final - t);
        auto result = ctrlpp::detail::rk45_step(exponential_decay{}, x, u, dt_step, cfg);
        x = result.x_next;
        t += result.dt_actual;
        dt = result.dt_next;
    }

    double exact = std::exp(-1.0);
    CHECK_THAT(x(0), Catch::Matchers::WithinAbs(exact, 1e-6));
}
