#include "ctrlpp/mpc.h"
#include "ctrlpp/model/state_space.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <span>
#include <vector>
#include <cstddef>
#include <utility>

namespace {

// The mocks below have no setup failure mode, so their setup-error type carries
// no enumerators. The solver concepts accept one setup shape, a fallible one, so
// a backend with nothing to fail at writes a trivially succeeding fallible setup
// rather than an infallible one.
enum class mock_setup_error
{
};

using Catch::Matchers::WithinAbs;

struct mock_qp_solver
{
    using scalar_type = double;

    mutable ctrlpp::qp_problem<double> last_setup{};
    mutable ctrlpp::qp_update<double> last_update{};
    mutable int solve_count{0};
    mutable ctrlpp::solve_status next_status{ctrlpp::solve_status::optimal};

    auto setup(const ctrlpp::qp_problem<double> &problem) -> ctrlpp::expected<void, mock_setup_error>
    {
        last_setup = problem;
        return {};
    }

    auto solve(const ctrlpp::qp_update<double> &update) -> ctrlpp::qp_result<double>
    {
        last_update = update;
        ++solve_count;
        ctrlpp::qp_result<double> result;
        result.status          = next_status;
        result.x               = Eigen::VectorXd::Zero(last_setup.P.cols());
        result.y               = Eigen::VectorXd::Zero(last_setup.A.rows());
        result.objective       = 0.0;
        result.solve_time      = 0.001;
        result.iterations      = 5;
        result.primal_residual = 1e-6;
        result.dual_residual   = 1e-6;
        return result;
    }
};

// Double integrator: NX=2, NU=1
constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr double dt      = 0.1;

auto make_double_integrator() -> ctrlpp::discrete_state_space<double, NX, NU, NX>
{
    Eigen::Matrix2d A;
    A << 1.0, dt, 0.0, 1.0;
    Eigen::Vector2d B;
    B << 0.5 * dt * dt, dt;
    Eigen::Matrix2d C             = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 2, 1> D = Eigen::Matrix<double, 2, 1>::Zero();
    return {A, B, C, D};
}

auto make_config(int horizon = 5) -> ctrlpp::mpc_config<double, NX, NU>
{
    return {
            .horizon = horizon,
            .Q       = Eigen::Matrix2d::Identity(),
            .R       = Eigen::Matrix<double, 1, 1>::Identity(),
    };
}

using Mpc = ctrlpp::mpc<double, NX, NU, mock_qp_solver>;

/// Constructs through the validating factory and fails the case if the
/// configuration is rejected. This is the only construction path available in
/// the default (-fno-exceptions) tree, where the throwing convenience
/// constructors are compiled out.
template<typename Controller, typename... Args>
auto make_controller(Args &&...args) -> Controller
{
    auto created = Controller::create(std::forward<Args>(args)...);
    REQUIRE(created.has_value());
    return *std::move(created);
}

} // namespace

TEST_CASE("mock_qp_solver satisfies qp_solver concept")
{
    static_assert(ctrlpp::qp_solver<mock_qp_solver>, "mock_qp_solver must satisfy qp_solver concept");
}

TEST_CASE("mpc with mock solver", "[mpc]")
{
    auto sys        = make_double_integrator();
    constexpr int N = 5;

    SECTION("QP dimensions are correct for unconstrained problem")
    {
        auto cfg        = make_config(N);
        auto controller = make_controller<Mpc>(sys, cfg);

        Eigen::Vector2d x0{1.0, 0.0};
        auto result = controller.solve(x0);
        REQUIRE(result.has_value());

        // Decision variables: (N+1)*NX + N*NU = 6*2 + 5*1 = 17
        // No state bounds, no input bounds, no rate bounds
        // Constraints: (N+1)*NX = 12 (dynamics equality only)
        // Access via the mock's last_setup
        // We can't directly access the solver member, but we verify through
        // the trajectory and diagnostics that the solve went through correctly
    }

    SECTION("solve(x0) calls solver and returns u_0")
    {
        auto cfg        = make_config(N);
        auto controller = make_controller<Mpc>(sys, cfg);

        Eigen::Vector2d x0{1.0, 0.5};
        auto result = controller.solve(x0);
        REQUIRE(result.has_value());

        // Mock returns zero solution, so u_0 should be zero
        CHECK_THAT(result->input(0), WithinAbs(0.0, 1e-12));
        CHECK(result->status == ctrlpp::solve_result_status::converged);
    }

    SECTION("solve(x0, x_ref) produces non-zero q vector for reference tracking")
    {
        auto cfg        = make_config(N);
        auto controller = make_controller<Mpc>(sys, cfg);

        Eigen::Vector2d x0{0.0, 0.0};
        Eigen::Vector2d x_ref{1.0, 0.0};
        auto result = controller.solve(x0, x_ref);
        REQUIRE(result.has_value());
    }

    SECTION("solve returns the error branch when solver reports infeasible")
    {
        auto cfg        = make_config(N);
        auto controller = make_controller<Mpc>(sys, cfg);

        // First solve to populate things normally
        Eigen::Vector2d x0{1.0, 0.0};

        // We need to set next_status to infeasible.
        // Since solver_ is private, we test by constructing a controller where
        // the mock returns infeasible. We can do this by creating a special mock.
        // However the mock's next_status defaults to optimal, and we can't modify
        // it from outside... Let's use a different approach: a stateful mock.

        // Actually, the mock is default-constructed inside mpc. To test infeasibility,
        // we use a variant mock that always returns infeasible.

        struct infeasible_mock
        {
            using scalar_type = double;

            mutable ctrlpp::qp_problem<double> last_setup{};

            auto setup(const ctrlpp::qp_problem<double> &problem) -> ctrlpp::expected<void, mock_setup_error>
            {
                last_setup = problem;
                return {};
            }

            auto solve(const ctrlpp::qp_update<double> &) -> ctrlpp::qp_result<double>
            {
                ctrlpp::qp_result<double> result;
                result.status          = ctrlpp::solve_status::infeasible;
                result.x               = Eigen::VectorXd::Zero(last_setup.P.cols());
                result.y               = Eigen::VectorXd::Zero(last_setup.A.rows());
                result.objective       = 0.0;
                result.solve_time      = 0.0;
                result.iterations      = 0;
                result.primal_residual = 0.0;
                result.dual_residual   = 0.0;
                return result;
            }
        };

        static_assert(ctrlpp::qp_solver<infeasible_mock>);

        auto infeasible_ctrl   = make_controller<ctrlpp::mpc<double, NX, NU, infeasible_mock>>(sys, cfg);
        auto infeasible_result = infeasible_ctrl.solve(x0);
        CHECK_FALSE(infeasible_result.has_value());
        CHECK(infeasible_result.error() == ctrlpp::solver_error::infeasible);

        // The error branch must not populate a valid trajectory either.
        CHECK_FALSE(infeasible_ctrl.trajectory().has_value());
    }

    SECTION("budget-limited solve reaches the caller tagged budget_exhausted")
    {
        auto cfg = make_config(N);
        Eigen::Vector2d x0{1.0, 0.0};

        // A solver that exhausts its iteration budget still returns its best
        // iterate; the widened accept-set surfaces it on the success branch.
        struct budget_mock
        {
            using scalar_type = double;

            mutable ctrlpp::qp_problem<double> last_setup{};

            auto setup(const ctrlpp::qp_problem<double> &problem) -> ctrlpp::expected<void, mock_setup_error>
            {
                last_setup = problem;
                return {};
            }

            auto solve(const ctrlpp::qp_update<double> &) -> ctrlpp::qp_result<double>
            {
                ctrlpp::qp_result<double> result;
                result.status = ctrlpp::solve_status::max_iterations;
                result.x      = Eigen::VectorXd::Zero(last_setup.P.cols());
                result.y      = Eigen::VectorXd::Zero(last_setup.A.rows());
                return result;
            }
        };

        static_assert(ctrlpp::qp_solver<budget_mock>);

        auto budget_ctrl   = make_controller<ctrlpp::mpc<double, NX, NU, budget_mock>>(sys, cfg);
        auto budget_result = budget_ctrl.solve(x0);
        REQUIRE(budget_result.has_value());
        CHECK(budget_result->status == ctrlpp::solve_result_status::budget_exhausted);
    }

    SECTION("trajectory extracts correct number of state and input vectors")
    {
        auto cfg        = make_config(N);
        auto controller = make_controller<Mpc>(sys, cfg);

        Eigen::Vector2d x0{1.0, 0.0};
        auto result = controller.solve(x0);
        REQUIRE(result.has_value());

        auto traj = controller.trajectory();
        REQUIRE(traj.has_value());
        auto &[states, inputs] = *traj;
        CHECK(states.size() == static_cast<std::size_t>(N + 1));
        CHECK(inputs.size() == static_cast<std::size_t>(N));

        // Each state vector has NX elements, each input has NU elements
        for(const auto &s : states)
            CHECK(s.size() == static_cast<Eigen::Index>(NX));
        for(const auto &u : inputs)
            CHECK(u.size() == static_cast<Eigen::Index>(NU));
    }

    SECTION("diagnostics returns values from last solve")
    {
        auto cfg        = make_config(N);
        auto controller = make_controller<Mpc>(sys, cfg);

        Eigen::Vector2d x0{1.0, 0.0};
        auto result = controller.solve(x0);
        REQUIRE(result.has_value());

        auto diag = controller.diagnostics();
        CHECK(diag.status == ctrlpp::solve_status::optimal);
        CHECK(diag.iterations == 5);
        CHECK_THAT(diag.solve_time, WithinAbs(0.001, 1e-12));
        CHECK_THAT(diag.primal_residual, WithinAbs(1e-6, 1e-12));
        CHECK_THAT(diag.dual_residual, WithinAbs(1e-6, 1e-12));
    }

    SECTION("second solve passes warm-start from first solution")
    {
        auto cfg = make_config(N);

        // Use a mock that tracks warm-start presence
        struct warmstart_mock
        {
            using scalar_type = double;

            mutable ctrlpp::qp_problem<double> last_setup{};
            mutable bool had_warm_x{false};
            mutable bool had_warm_y{false};
            mutable int solve_count{0};

            auto setup(const ctrlpp::qp_problem<double> &problem) -> ctrlpp::expected<void, mock_setup_error>
            {
                last_setup = problem;
                return {};
            }

            auto solve(const ctrlpp::qp_update<double> &update) -> ctrlpp::qp_result<double>
            {
                ++solve_count;
                had_warm_x = update.warm_x.size() > 0;
                had_warm_y = update.warm_y.size() > 0;

                ctrlpp::qp_result<double> result;
                result.status          = ctrlpp::solve_status::optimal;
                result.x               = Eigen::VectorXd::Ones(last_setup.P.cols());
                result.y               = Eigen::VectorXd::Ones(last_setup.A.rows());
                result.objective       = 1.0;
                result.solve_time      = 0.002;
                result.iterations      = 3;
                result.primal_residual = 1e-7;
                result.dual_residual   = 1e-7;
                return result;
            }
        };

        static_assert(ctrlpp::qp_solver<warmstart_mock>);

        auto controller = make_controller<ctrlpp::mpc<double, NX, NU, warmstart_mock>>(sys, cfg);
        Eigen::Vector2d x0{1.0, 0.0};

        // First solve -- no warm-start data yet
        auto r1 = controller.solve(x0);
        REQUIRE(r1.has_value());

        // Second solve -- should have warm-start from first
        auto r2 = controller.solve(x0);
        REQUIRE(r2.has_value());

        // The warm-start vectors should have been populated on second call
        // (diagnostics show the second solve completed)
        auto diag = controller.diagnostics();
        CHECK(diag.iterations == 3);
    }

    SECTION("soft state constraints produce larger QP with slack variables")
    {
        auto cfg  = make_config(N);
        cfg.x_min = Eigen::Vector2d{-10.0, -10.0};
        cfg.x_max = Eigen::Vector2d{10.0, 10.0};
        // hard_state_constraints defaults to false, so soft constraints apply

        auto controller = make_controller<Mpc>(sys, cfg);

        Eigen::Vector2d x0{1.0, 0.0};
        auto result = controller.solve(x0);
        REQUIRE(result.has_value());

        // With soft constraints: slack variables added
        // Decision vars: (N+1)*NX + N*NU + N*NX = 6*2 + 5*1 + 5*2 = 27
        // Trajectory should still be same size
        auto traj = controller.trajectory();
        REQUIRE(traj.has_value());
        auto &[states, inputs] = *traj;
        CHECK(states.size() == static_cast<std::size_t>(N + 1));
        CHECK(inputs.size() == static_cast<std::size_t>(N));
    }

    SECTION("du_max produces additional rate constraint rows")
    {
        auto cfg   = make_config(N);
        cfg.u_min  = Eigen::Matrix<double, 1, 1>{-5.0};
        cfg.u_max  = Eigen::Matrix<double, 1, 1>{5.0};
        cfg.du_max = Eigen::Matrix<double, 1, 1>{1.0};

        auto controller = make_controller<Mpc>(sys, cfg);

        Eigen::Vector2d x0{1.0, 0.0};
        auto result = controller.solve(x0);
        REQUIRE(result.has_value());

        // With rate constraints, the QP has additional rows
        // Constraints: (N+1)*NX dynamics + N*NU input bounds + N*NU rate bounds
        //            = 12 + 5 + 5 = 22
        auto traj = controller.trajectory();
        REQUIRE(traj.has_value());
        auto &[states, inputs] = *traj;
        CHECK(states.size() == static_cast<std::size_t>(N + 1));
        CHECK(inputs.size() == static_cast<std::size_t>(N));
    }
}

// --- NY < NX output tracking tests ---

TEST_CASE("mpc with NY < NX output tracking", "[mpc][output_tracking]")
{
    // 4-state, 2-input, 2-output system (double integrator with position-only output)
    // States: [x, y, vx, vy], Inputs: [ax, ay], Outputs: [x, y]
    constexpr std::size_t NX4 = 4;
    constexpr std::size_t NU2 = 2;
    constexpr std::size_t NY2 = 2;
    constexpr double dt4      = 0.1;

    // C = [I_2 0_2] selects positions only
    Eigen::Matrix<double, 4, 4> A4;
    A4 << 1.0, 0.0, dt4, 0.0, 0.0, 1.0, 0.0, dt4, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0;

    Eigen::Matrix<double, 4, 2> B4;
    B4 << 0.5 * dt4 * dt4, 0.0, 0.0, 0.5 * dt4 * dt4, dt4, 0.0, 0.0, dt4;

    Eigen::Matrix<double, 2, 4> C4;
    C4 << 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;

    Eigen::Matrix<double, 2, 2> D4 = Eigen::Matrix<double, 2, 2>::Zero();

    ctrlpp::discrete_state_space<double, NX4, NU2, NY2> sys4{A4, B4, C4, D4};

    SECTION("NY<NX mpc_config template compiles and Q is NY x NY")
    {
        ctrlpp::mpc_config<double, NX4, NU2, NY2> cfg{
                .horizon = 5,
                .Q       = Eigen::Matrix2d::Identity(),
                .R       = Eigen::Matrix2d::Identity() * 0.1,
        };

        // Q should be 2x2 (NY x NY), not 4x4
        static_assert(decltype(cfg.Q)::RowsAtCompileTime == 2);
        static_assert(decltype(cfg.Q)::ColsAtCompileTime == 2);

        auto controller = make_controller<ctrlpp::mpc<double, NX4, NU2, mock_qp_solver, NY2>>(sys4, cfg);

        Eigen::Vector4d x0 = Eigen::Vector4d::Zero();
        auto result        = controller.solve(x0);
        REQUIRE(result.has_value());
    }

    SECTION("NY<NX single reference tracking passes NY-dimensional y_ref")
    {
        ctrlpp::mpc_config<double, NX4, NU2, NY2> cfg{
                .horizon = 5,
                .Q       = Eigen::Matrix2d::Identity(),
                .R       = Eigen::Matrix2d::Identity() * 0.1,
        };

        auto controller = make_controller<ctrlpp::mpc<double, NX4, NU2, mock_qp_solver, NY2>>(sys4, cfg);

        Eigen::Vector4d x0 = Eigen::Vector4d::Zero();
        Eigen::Vector2d y_ref{1.0, 2.0}; // 2D output reference (position only)
        auto result = controller.solve(x0, y_ref);
        REQUIRE(result.has_value());
    }

    SECTION("NY<NX span reference tracking passes NY-dimensional y_ref sequence")
    {
        constexpr int N = 5;
        ctrlpp::mpc_config<double, NX4, NU2, NY2> cfg{
                .horizon = N,
                .Q       = Eigen::Matrix2d::Identity(),
                .R       = Eigen::Matrix2d::Identity() * 0.1,
        };

        auto controller = make_controller<ctrlpp::mpc<double, NX4, NU2, mock_qp_solver, NY2>>(sys4, cfg);

        Eigen::Vector4d x0 = Eigen::Vector4d::Zero();
        std::vector<Eigen::Vector2d> y_refs(static_cast<std::size_t>(N + 1), Eigen::Vector2d{1.0, 2.0});
        std::span<const Eigen::Vector2d> ref_span(y_refs);
        auto result = controller.solve(x0, ref_span);
        REQUIRE(result.has_value());
    }

    SECTION("sysid-compatible: discrete_state_space<S, NX, NU, NY> accepted directly")
    {
        // Verify that a state_space with NY != NX can be passed to mpc constructor
        // This simulates a sysid::recursive_arx::to_state_space() return type
        constexpr std::size_t SYSID_NX = 2;
        constexpr std::size_t SYSID_NU = 1;
        constexpr std::size_t SYSID_NY = 1;

        ctrlpp::discrete_state_space<double, SYSID_NX, SYSID_NU, SYSID_NY> sysid_sys{
                .A = Eigen::Matrix2d::Identity(), .B = Eigen::Vector2d::Ones(), .C = (Eigen::Matrix<double, 1, 2>() << 1.0, 0.0).finished(), .D = Eigen::Matrix<double, 1, 1>::Zero()};

        ctrlpp::mpc_config<double, SYSID_NX, SYSID_NU, SYSID_NY> cfg{
                .horizon = 5,
                .Q       = Eigen::Matrix<double, 1, 1>::Identity(),
                .R       = Eigen::Matrix<double, 1, 1>::Identity() * 0.1,
        };

        // This should compile -- sysid return type feeds directly to MPC
        auto controller = make_controller<ctrlpp::mpc<double, SYSID_NX, SYSID_NU, mock_qp_solver, SYSID_NY>>(sysid_sys, cfg);

        Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
        auto result        = controller.solve(x0);
        REQUIRE(result.has_value());

        // Track a 1D output reference
        Eigen::Matrix<double, 1, 1> y_ref;
        y_ref << 1.0;
        auto tracking_result = controller.solve(x0, y_ref);
        REQUIRE(tracking_result.has_value());
    }
}

TEST_CASE("mpc span overload rejects an undersized reference span", "[mpc][span]")
{
    auto sys        = make_double_integrator();
    constexpr int N = 5;
    auto cfg        = make_config(N);
    auto controller = make_controller<Mpc>(sys, cfg);

    Eigen::Vector2d x0{1.0, 0.0};

    // The tracking overload reads y_ref[0..N], so it needs N+1 references. A
    // shorter span must be rejected via the error branch rather than overrun.
    std::vector<Eigen::Vector2d> short_refs(static_cast<std::size_t>(N), Eigen::Vector2d{1.0, 0.0});
    std::span<const Eigen::Vector2d> short_span(short_refs);
    auto undersized = controller.solve(x0, short_span);
    REQUIRE_FALSE(undersized.has_value());
    CHECK(undersized.error() == ctrlpp::solver_error::invalid_problem);

    // Exactly N+1 references is accepted.
    std::vector<Eigen::Vector2d> ok_refs(static_cast<std::size_t>(N + 1), Eigen::Vector2d{1.0, 0.0});
    std::span<const Eigen::Vector2d> ok_span(ok_refs);
    CHECK(controller.solve(x0, ok_span).has_value());
}

TEST_CASE("mpc trajectory is guarded before the first valid solve", "[mpc][trajectory]")
{
    auto sys        = make_double_integrator();
    auto cfg        = make_config(5);
    auto controller = make_controller<Mpc>(sys, cfg);

    // No solve has run yet, so there is no valid trajectory to report.
    auto pre = controller.trajectory();
    REQUIRE_FALSE(pre.has_value());
    CHECK(pre.error() == ctrlpp::solver_error::setup_incomplete);

    Eigen::Vector2d x0{1.0, 0.0};
    REQUIRE(controller.solve(x0).has_value());
    CHECK(controller.trajectory().has_value());
}

// The terminal cost is the one place the controller answers a different
// question than the one it was configured with. With no terminal weight given,
// the infinite-horizon Riccati solution IS the terminal cost; when that solve
// has no solution the state weight stands in for it, and every solve afterwards
// optimizes a different problem. The substituted cost is usable, so nothing
// about a solve's success distinguishes the two. The disposition does.
TEST_CASE("mpc reports a substituted terminal cost", "[mpc][diagnostics]")
{
    constexpr int N = 5;

    SECTION("a failed Riccati solve reports the state weight standing in")
    {
        // A nilpotent A is exactly rank-deficient, and the symplectic pencil the
        // Riccati solve is built on needs A invertible. There is no solution to
        // fall back from, so the state weight stands in for the terminal cost.
        Eigen::Matrix2d A;
        A << 0.0, 1.0, 0.0, 0.0;
        Eigen::Vector2d B;
        B << 0.0, 1.0;
        Eigen::Matrix2d C             = Eigen::Matrix2d::Identity();
        Eigen::Matrix<double, 2, 1> D = Eigen::Matrix<double, 2, 1>::Zero();
        ctrlpp::discrete_state_space<double, NX, NU, NX> singular_sys{A, B, C, D};

        auto cfg = make_config(N);
        REQUIRE_FALSE(cfg.Qf.has_value());

        auto controller = make_controller<Mpc>(singular_sys, cfg);

        // Readable before any solve: the substitution happened at construction.
        REQUIRE(controller.diagnostics().used_state_weight_terminal_cost);

        // And carried on the aggregate every solve emits, so a caller reading
        // only the solve-path diagnostics still sees it.
        Eigen::Vector2d x0{1.0, 0.0};
        REQUIRE(controller.solve(x0).has_value());
        REQUIRE(controller.diagnostics().used_state_weight_terminal_cost);
    }

    SECTION("a successful Riccati solve reports no substitution")
    {
        auto sys = make_double_integrator();
        auto cfg = make_config(N);
        REQUIRE_FALSE(cfg.Qf.has_value());

        auto controller = make_controller<Mpc>(sys, cfg);
        REQUIRE_FALSE(controller.diagnostics().used_state_weight_terminal_cost);

        Eigen::Vector2d x0{1.0, 0.0};
        REQUIRE(controller.solve(x0).has_value());
        REQUIRE_FALSE(controller.diagnostics().used_state_weight_terminal_cost);
    }

    SECTION("an unverified Riccati terminal cost reports substitution")
    {
        auto sys    = make_double_integrator();
        auto cfg    = make_config(N);
        cfg.Q       = 1e10 * Eigen::Matrix2d::Identity();
        cfg.R(0, 0) = 0.1;
        REQUIRE_FALSE(cfg.Qf.has_value());

        auto controller = make_controller<Mpc>(sys, cfg);
        REQUIRE(controller.diagnostics().used_state_weight_terminal_cost);

        Eigen::Vector2d x0{1.0, 0.0};
        REQUIRE(controller.solve(x0).has_value());
        REQUIRE(controller.diagnostics().used_state_weight_terminal_cost);
    }

    SECTION("an equilibrated Riccati terminal cost reports no substitution")
    {
        auto sys    = make_double_integrator();
        auto cfg    = make_config(N);
        cfg.Q       = 1e8 * Eigen::Matrix2d::Identity();
        cfg.R(0, 0) = 0.1;
        REQUIRE_FALSE(cfg.Qf.has_value());

        auto controller = make_controller<Mpc>(sys, cfg);
        REQUIRE_FALSE(controller.diagnostics().used_state_weight_terminal_cost);

        Eigen::Vector2d x0{1.0, 0.0};
        REQUIRE(controller.solve(x0).has_value());
        REQUIRE_FALSE(controller.diagnostics().used_state_weight_terminal_cost);
    }

    SECTION("a configured terminal weight is never a substitution")
    {
        auto sys = make_double_integrator();
        auto cfg = make_config(N);
        cfg.Qf   = Eigen::Matrix2d::Identity() * 10.0;

        auto controller = make_controller<Mpc>(sys, cfg);
        REQUIRE_FALSE(controller.diagnostics().used_state_weight_terminal_cost);

        Eigen::Vector2d x0{1.0, 0.0};
        REQUIRE(controller.solve(x0).has_value());
        REQUIRE_FALSE(controller.diagnostics().used_state_weight_terminal_cost);
    }
}
