#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/state_space.h"

#include <catch2/catch_test_macros.hpp>
#include <rapidcheck/catch.h>
#include <rapidcheck.h>

#include <cstddef>

namespace {

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;

using OsqpMpc = ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver>;

auto bounded_double(double lo, double hi) -> rc::Gen<double>
{
    return rc::gen::map(rc::gen::inRange(0, 1000000), [lo, hi](int x) {
        return lo + (hi - lo) * (static_cast<double>(x) / 1000000.0);
    });
}

auto make_double_integrator(double dt)
    -> ctrlpp::discrete_state_space<double, NX, NU, NX>
{
    Eigen::Matrix2d A;
    A << 1.0, dt,
         0.0, 1.0;
    Eigen::Vector2d B;
    B << 0.5 * dt * dt,
         dt;
    Eigen::Matrix2d C = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 2, 1> D = Eigen::Matrix<double, 2, 1>::Zero();
    return {A, B, C, D};
}

}

TEST_CASE("mpc input constraints satisfied", "[mpc][property]")
{
    rc::prop("optimal input respects u_min/u_max bounds", [] {
        auto dt = *bounded_double(0.01, 0.1);
        auto sys = make_double_integrator(dt);

        auto u_lo = *bounded_double(-20.0, -1.0);
        auto u_hi = *bounded_double(1.0, 20.0);
        auto horizon = *rc::gen::inRange(3, 16);
        auto x0_0 = *bounded_double(-5.0, 5.0);
        auto x0_1 = *bounded_double(-5.0, 5.0);

        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = horizon,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
            .u_min = (Eigen::Matrix<double, 1, 1>() << u_lo).finished(),
            .u_max = (Eigen::Matrix<double, 1, 1>() << u_hi).finished(),
        };

        OsqpMpc controller(sys, cfg);
        Eigen::Vector2d x0{x0_0, x0_1};
        auto u = controller.solve(x0);

        if (u.has_value()) {
            RC_ASSERT((*u)(0) >= u_lo - 1e-6);
            RC_ASSERT((*u)(0) <= u_hi + 1e-6);

            auto [states, inputs] = controller.trajectory();
            for (const auto& ui : inputs) {
                RC_ASSERT(ui(0) >= u_lo - 1e-6);
                RC_ASSERT(ui(0) <= u_hi + 1e-6);
            }
        }
    });
}

TEST_CASE("mpc predicted states within bounds", "[mpc][property]")
{
    rc::prop("predicted states respect x_min/x_max with slack tolerance", [] {
        auto dt = *bounded_double(0.01, 0.1);
        auto sys = make_double_integrator(dt);

        auto x0_0 = *bounded_double(-2.0, 2.0);
        auto x0_1 = *bounded_double(-2.0, 2.0);

        // Ensure feasible: x_min < x0 < x_max
        auto margin = *bounded_double(1.0, 5.0);
        Eigen::Vector2d x_lo{x0_0 - margin, x0_1 - margin};
        Eigen::Vector2d x_hi{x0_0 + margin, x0_1 + margin};

        auto horizon = *rc::gen::inRange(3, 16);

        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = horizon,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
            .x_min = x_lo,
            .x_max = x_hi,
        };

        OsqpMpc controller(sys, cfg);
        Eigen::Vector2d x0{x0_0, x0_1};
        auto u = controller.solve(x0);

        if (u.has_value()) {
            auto [states, inputs] = controller.trajectory();
            // Soft constraints: allow slack up to penalty tolerance
            constexpr double slack_tol = 0.1;
            for (const auto& s : states) {
                RC_ASSERT(s(0) >= x_lo(0) - slack_tol);
                RC_ASSERT(s(0) <= x_hi(0) + slack_tol);
                RC_ASSERT(s(1) >= x_lo(1) - slack_tol);
                RC_ASSERT(s(1) <= x_hi(1) + slack_tol);
            }
        }
    });
}

TEST_CASE("mpc cost is non-negative", "[mpc][property]")
{
    rc::prop("cost from diagnostics is non-negative for PSD Q, R", [] {
        auto dt = *bounded_double(0.01, 0.1);
        auto sys = make_double_integrator(dt);

        // Build PSD Q via L^T * L
        auto l00 = *bounded_double(0.1, 5.0);
        auto l10 = *bounded_double(-2.0, 2.0);
        auto l11 = *bounded_double(0.1, 5.0);
        Eigen::Matrix2d L;
        L << l00, 0.0,
             l10, l11;
        Eigen::Matrix2d Q = L.transpose() * L;

        // PD R
        auto r_val = *bounded_double(0.01, 10.0);
        auto R = (Eigen::Matrix<double, 1, 1>() << r_val).finished();

        auto horizon = *rc::gen::inRange(3, 16);
        auto x0_0 = *bounded_double(-5.0, 5.0);
        auto x0_1 = *bounded_double(-5.0, 5.0);

        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = horizon,
            .Q = Q,
            .R = R,
        };

        OsqpMpc controller(sys, cfg);
        Eigen::Vector2d x0{x0_0, x0_1};
        auto u = controller.solve(x0);

        if (u.has_value()) {
            auto diag = controller.diagnostics();
            RC_ASSERT(diag.cost >= -1e-10);
        }
    });
}

TEST_CASE("mpc robustness - infeasible returns nullopt", "[mpc][property]")
{
    rc::prop("contradictory constraints return nullopt or slack solution, no crash", [] {
        auto dt = *bounded_double(0.01, 0.1);
        auto sys = make_double_integrator(dt);

        auto horizon = *rc::gen::inRange(3, 10);

        // Very tight state bounds around origin with large initial state
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = horizon,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
            .u_min = (Eigen::Matrix<double, 1, 1>() << -0.01).finished(),
            .u_max = (Eigen::Matrix<double, 1, 1>() << 0.01).finished(),
            .x_min = (Eigen::Vector2d() << -0.01, -0.01).finished(),
            .x_max = (Eigen::Vector2d() << 0.01, 0.01).finished(),
            .hard_state_constraints = true,
        };

        OsqpMpc controller(sys, cfg);
        Eigen::Vector2d x0{10.0, 10.0};

        // Must not crash -- either returns nullopt or a solution with slack
        auto result = controller.solve(x0);
        // With hard constraints and infeasible config, expect nullopt
        RC_ASSERT(!result.has_value());
    });
}
