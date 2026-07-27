#include "ctrlpp/mpc.h"
#include "ctrlpp/expected.h"
#include "ctrlpp/mpc/osqp_solver.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <type_traits>

namespace
{

using Catch::Matchers::WithinAbs;

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr double dt = 0.1;

auto make_double_integrator() -> ctrlpp::discrete_state_space<double, NX, NU, NX>
{
    Eigen::Matrix2d A;
    A << 1.0, dt, 0.0, 1.0;
    Eigen::Vector2d B;
    B << 0.5 * dt * dt, dt;
    Eigen::Matrix2d C = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 2, 1> D = Eigen::Matrix<double, 2, 1>::Zero();
    return {A, B, C, D};
}

using OsqpMpc = ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver>;

auto make_simple_qp() -> ctrlpp::qp_problem<double>
{
    Eigen::SparseMatrix<double> P(2, 2);
    P.insert(0, 0) = 1.0;
    P.insert(1, 1) = 1.0;
    P.makeCompressed();

    Eigen::SparseMatrix<double> A(2, 2);
    A.insert(0, 0) = 1.0;
    A.insert(1, 1) = 1.0;
    A.makeCompressed();

    return {.P = P, .q = Eigen::Vector2d::Zero(), .A = A, .l = Eigen::Vector2d::Constant(-1.0), .u = Eigen::Vector2d::Constant(1.0)};
}

} // namespace

TEST_CASE("osqp_solver try_setup reports setup failure as an expected error", "[mpc][osqp]")
{
    static_assert(std::is_same_v<decltype(std::declval<ctrlpp::osqp_solver&>().try_setup(std::declval<const ctrlpp::qp_problem<double>&>())), ctrlpp::expected<void, ctrlpp::osqp_setup_error>>,
                  "try_setup must return ctrlpp::expected<void, osqp_setup_error>");

    auto problem = make_simple_qp();

    SECTION("well-formed problem sets up successfully")
    {
        ctrlpp::osqp_solver solver;
        auto result = solver.try_setup(problem);
        REQUIRE(result.has_value());
    }

    SECTION("rejected solver settings surface as setup_failed instead of a throw")
    {
        // OSQP validates settings inside osqp_setup; a negative absolute
        // tolerance fails that validation deterministically.
        ctrlpp::osqp_solver solver{-1.0};
        auto result = solver.try_setup(problem);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == ctrlpp::osqp_setup_error::setup_failed);
    }

    SECTION("a setup that fails after OSQP has allocated leaves the solver reusable")
    {
        // osqp_setup publishes its solver pointer before the allocations that
        // can fail, and its error paths free nothing, so the caller owns the
        // partially built solver. A NaN cost weight passes OSQP's up-front data
        // validation and fails later, once that allocation has happened --
        // unlike the rejected-settings case above, which fails before it.
        // Releasing it is what keeps the solver usable for a second setup;
        // under ASan/LSan, dropping the pointer instead shows up as a leak.
        auto nan_problem = make_simple_qp();
        nan_problem.P.coeffRef(0, 0) = std::numeric_limits<double>::quiet_NaN();
        nan_problem.P.makeCompressed();

        ctrlpp::osqp_solver solver;
        REQUIRE_FALSE(solver.try_setup(nan_problem).has_value());

        REQUIRE(solver.try_setup(make_simple_qp()).has_value());

        ctrlpp::qp_update<double> update{.q = Eigen::Vector2d::Zero(), .l = Eigen::Vector2d::Constant(-1.0), .u = Eigen::Vector2d::Constant(1.0), .warm_x = {}, .warm_y = {}};
        CHECK(solver.solve(update).status == ctrlpp::solve_status::optimal);
    }

    SECTION("solve_status stays a plain value on the solve path")
    {
        static_assert(std::is_same_v<decltype(ctrlpp::qp_result<double>{}.status), ctrlpp::solve_status>, "qp_result::status must stay a plain solve_status value");

        ctrlpp::osqp_solver solver;
        REQUIRE(solver.try_setup(problem).has_value());

        ctrlpp::qp_update<double> update{.q = Eigen::Vector2d::Zero(), .l = Eigen::Vector2d::Constant(-1.0), .u = Eigen::Vector2d::Constant(1.0), .warm_x = {}, .warm_y = {}};
        auto result = solver.solve(update);
        CHECK(result.status == ctrlpp::solve_status::optimal);
    }
}

TEST_CASE("mpc accepts a preset-injected solver", "[mpc][osqp]")
{
    auto sys = make_double_integrator();
    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = 10,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
    };

    // Third constructor argument injects a pre-configured solver; qp_preset::speed
    // builds it with polishing off. The closed-loop regulation must still hold.
    auto controller_result = OsqpMpc::create(sys, cfg, ctrlpp::osqp_solver{ctrlpp::qp_preset::speed});
    REQUIRE(controller_result.has_value());
    auto& controller = *controller_result;

    Eigen::Vector2d x{1.0, 0.0};
    for(int step = 0; step < 50; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = sys.A * x + sys.B * u.value().input;
    }
    CHECK(x.norm() < 0.1);
}

TEST_CASE("mpc with OSQP solver", "[mpc][osqp]")
{
    auto sys = make_double_integrator();
    constexpr int N = 10;

    SECTION("regulation drives state toward origin")
    {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        };

        auto controller_result = OsqpMpc::create(sys, cfg);
        REQUIRE(controller_result.has_value());
        auto& controller = *controller_result;

        Eigen::Vector2d x{1.0, 0.0};
        double prev_norm = x.norm();

        // Simulate closed-loop for 50 steps
        int monotonic_after = 5; // allow transient in first few steps
        for(int step = 0; step < 50; ++step)
        {
            auto u = controller.solve(x);
            REQUIRE(u.has_value());
            x = sys.A * x + sys.B * u.value().input;

            double norm = x.norm();
            if(step >= monotonic_after)
            {
                CHECK(norm <= prev_norm + 1e-6);
            }
            prev_norm = norm;
        }

        CHECK(x.norm() < 0.1);
    }

    SECTION("setpoint tracking converges to reference")
    {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        };

        auto controller_result = OsqpMpc::create(sys, cfg);
        REQUIRE(controller_result.has_value());
        auto& controller = *controller_result;

        Eigen::Vector2d x{0.0, 0.0};
        Eigen::Vector2d x_ref{2.0, 0.0};

        // First solve should produce non-zero control pushing toward reference
        auto u0 = controller.solve(x, x_ref);
        REQUIRE(u0.has_value());
        CHECK(std::abs(u0->input(0)) > 1e-6);

        // Simulate closed-loop for 50 steps
        for(int step = 0; step < 50; ++step)
        {
            auto u = controller.solve(x, x_ref);
            REQUIRE(u.has_value());
            x = sys.A * x + sys.B * u.value().input;
        }

        CHECK((x - x_ref).norm() < 0.1);
    }

    SECTION("input bounds are respected")
    {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
            .u_min = (Eigen::Matrix<double, 1, 1>() << -0.5).finished(),
            .u_max = (Eigen::Matrix<double, 1, 1>() << 0.5).finished(),
        };

        auto controller_result = OsqpMpc::create(sys, cfg);
        REQUIRE(controller_result.has_value());
        auto& controller = *controller_result;

        Eigen::Vector2d x{5.0, 0.0}; // large initial state to push solver hard

        for(int step = 0; step < 20; ++step)
        {
            auto u = controller.solve(x);
            REQUIRE(u.has_value());
            CHECK(u->input(0) >= -0.5 - 1e-4);
            CHECK(u->input(0) <= 0.5 + 1e-4);
            x = sys.A * x + sys.B * u.value().input;
        }
    }

    SECTION("soft state constraints allow feasible solution from outside bounds")
    {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
            .x_min = (Eigen::Vector2d() << -5.0, -5.0).finished(),
            .x_max = (Eigen::Vector2d() << 5.0, 5.0).finished(),
            .hard_state_constraints = false,
        };

        Eigen::Vector2d x0{10.0, 0.0}; // outside state bounds

        auto controller_result = OsqpMpc::create(sys, cfg);
        REQUIRE(controller_result.has_value());
        auto& controller = *controller_result;
        auto u = controller.solve(x0);
        REQUIRE(u.has_value());
    }

    SECTION("hard state constraints make out-of-bounds initial state infeasible")
    {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
            .x_min = (Eigen::Vector2d() << -5.0, -5.0).finished(),
            .x_max = (Eigen::Vector2d() << 5.0, 5.0).finished(),
            .hard_state_constraints = true,
        };

        Eigen::Vector2d x0{10.0, 0.0}; // violates hard state bounds

        auto controller_result = OsqpMpc::create(sys, cfg);
        REQUIRE(controller_result.has_value());
        auto& controller = *controller_result;
        auto result = controller.solve(x0);
        CHECK_FALSE(result.has_value());
    }

    SECTION("warm-starting reduces iterations on second solve")
    {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        };

        auto controller_result = OsqpMpc::create(sys, cfg);
        REQUIRE(controller_result.has_value());
        auto& controller = *controller_result;

        Eigen::Vector2d x0{1.0, 0.0};

        // First solve (cold)
        auto r1 = controller.solve(x0);
        REQUIRE(r1.has_value());
        int iter1 = controller.diagnostics().iterations;

        // Simulate one step forward so state changes slightly
        Eigen::Vector2d x1 = sys.A * x0 + sys.B * r1.value().input;

        // Second solve with nearby state (warm-started from first)
        auto r2 = controller.solve(x1);
        REQUIRE(r2.has_value());
        int iter2 = controller.diagnostics().iterations;

        CHECK(iter2 <= iter1);
    }

    SECTION("diagnostics populated after successful solve")
    {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        };

        auto controller_result = OsqpMpc::create(sys, cfg);
        REQUIRE(controller_result.has_value());
        auto& controller = *controller_result;

        Eigen::Vector2d x0{1.0, 0.0};
        auto result = controller.solve(x0);
        REQUIRE(result.has_value());

        auto diag = controller.diagnostics();
        CHECK(diag.iterations > 0);
        CHECK(diag.solve_time > 0.0);
        CHECK(diag.cost >= 0.0);
    }

    SECTION("trajectory returns consistent states and inputs")
    {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        };

        auto controller_result = OsqpMpc::create(sys, cfg);
        REQUIRE(controller_result.has_value());
        auto& controller = *controller_result;

        Eigen::Vector2d x0{1.0, 0.0};
        auto result = controller.solve(x0);
        REQUIRE(result.has_value());

        auto traj = controller.trajectory();
        REQUIRE(traj.has_value());
        auto& [states, inputs] = *traj;
        REQUIRE(states.size() == static_cast<std::size_t>(N + 1));
        REQUIRE(inputs.size() == static_cast<std::size_t>(N));

        // states[0] should match x0
        CHECK_THAT(states[0](0), WithinAbs(x0(0), 1e-3));
        CHECK_THAT(states[0](1), WithinAbs(x0(1), 1e-3));

        // States should propagate approximately via dynamics
        for(int k = 0; k < N; ++k)
        {
            Eigen::Vector2d x_next = sys.A * states[static_cast<std::size_t>(k)] + sys.B * inputs[static_cast<std::size_t>(k)];
            CHECK_THAT(states[static_cast<std::size_t>(k + 1)](0), WithinAbs(x_next(0), 1e-2));
            CHECK_THAT(states[static_cast<std::size_t>(k + 1)](1), WithinAbs(x_next(1), 1e-2));
        }
    }

    SECTION("rate constraints limit control change between solves")
    {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
            .u_min = (Eigen::Matrix<double, 1, 1>() << -5.0).finished(),
            .u_max = (Eigen::Matrix<double, 1, 1>() << 5.0).finished(),
            .du_max = (Eigen::Matrix<double, 1, 1>() << 0.2).finished(),
        };

        auto controller_result = OsqpMpc::create(sys, cfg);
        REQUIRE(controller_result.has_value());
        auto& controller = *controller_result;

        Eigen::Vector2d x{5.0, 0.0};
        double u_prev = 0.0; // initial u_prev_ is zero in mpc

        for(int step = 0; step < 15; ++step)
        {
            auto u = controller.solve(x);
            REQUIRE(u.has_value());
            double u_cur = u->input(0);

            CHECK(std::abs(u_cur - u_prev) <= 0.2 + 1e-2);

            u_prev = u_cur;
            x = sys.A * x + sys.B * u.value().input;
        }
    }

    SECTION("DARE default Qf solves without explicit Qf")
    {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
            // Qf not specified -- DARE should compute it internally
        };

        auto controller_result = OsqpMpc::create(sys, cfg);
        REQUIRE(controller_result.has_value());
        auto& controller = *controller_result;

        Eigen::Vector2d x0{1.0, 0.0};
        auto result = controller.solve(x0);
        REQUIRE(result.has_value());

        auto diag = controller.diagnostics();
        CHECK(diag.status == ctrlpp::solve_status::optimal);
    }

    SECTION("set_applied_input re-anchors the rate constraint")
    {
        ctrlpp::mpc_config<double, NX, NU> cfg{
            .horizon = N,
            .Q = Eigen::Matrix2d::Identity(),
            .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
            .u_min = (Eigen::Matrix<double, 1, 1>() << -5.0).finished(),
            .u_max = (Eigen::Matrix<double, 1, 1>() << 5.0).finished(),
            .du_max = (Eigen::Matrix<double, 1, 1>() << 0.2).finished(),
        };

        auto controller_result = OsqpMpc::create(sys, cfg);
        REQUIRE(controller_result.has_value());
        auto& controller = *controller_result;

        Eigen::Vector2d x{5.0, 0.0};

        // First solve anchors the rate window at the internal u_prev (zero).
        auto u0 = controller.solve(x);
        REQUIRE(u0.has_value());

        // Record a commanded input the caller actually applied; the next solve
        // must keep its input within du_max of that recorded value, not of the
        // solver's own last iterate.
        Eigen::Matrix<double, 1, 1> applied;
        applied << -3.0;
        controller.set_applied_input(applied);

        auto u1 = controller.solve(x);
        REQUIRE(u1.has_value());
        CHECK(std::abs(u1->input(0) - applied(0)) <= 0.2 + 1e-2);
    }
}
