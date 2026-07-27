#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/invariant_set.h"
#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/mpc/terminal_set.h"

#include "ctrlpp/control/dare.h"
#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <cstddef>

namespace
{

using Catch::Matchers::WithinAbs;

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr double dt = 0.1;

// Unbounded state box: the terminal ellipsoid is then constrained by the input faces alone.
auto no_state_min() -> Eigen::Vector2d
{
    return Eigen::Vector2d::Constant(-std::numeric_limits<double>::infinity());
}
auto no_state_max() -> Eigen::Vector2d
{
    return Eigen::Vector2d::Constant(std::numeric_limits<double>::infinity());
}

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

} // namespace

TEST_CASE("terminal_ingredients on stable double integrator", "[terminal_set]")
{
    auto sys = make_double_integrator();
    Eigen::Matrix2d Q = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 0.1;
    Eigen::Matrix<double, 1, 1> u_min;
    u_min << -1.0;
    Eigen::Matrix<double, 1, 1> u_max;
    u_max << 1.0;

    auto result = ctrlpp::terminal_ingredients<double, NX, NU>(sys.A, sys.B, Q, R, u_min, u_max, no_state_min(), no_state_max());

    REQUIRE(result.has_value());

    // Qf should be positive semi-definite
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(result->Qf);
    for(int i = 0; i < 2; ++i)
        CHECK(eig.eigenvalues()(i) >= -1e-10);

    // alpha should be positive and finite
    CHECK(result->set.alpha > 0.0);
    CHECK(std::isfinite(result->set.alpha));
}

TEST_CASE("compute_ellipsoidal_set with known DARE solution", "[terminal_set]")
{
    // Simple 2x2 system: P = I, K = [k1, k2]
    // alpha_i = min(u_max^2, u_min^2) / (k_i' P^{-1} k_i)
    // With P = I, P^{-1} = I, so alpha = u_bound^2 / ||k||^2
    Eigen::Matrix2d P = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 1, 2> K;
    K << 0.5, 0.5;
    Eigen::Matrix<double, 1, 1> u_min;
    u_min << -1.0;
    Eigen::Matrix<double, 1, 1> u_max;
    u_max << 1.0;

    auto eset = ctrlpp::compute_ellipsoidal_set<double, NX, NU>(P, K, u_min, u_max, no_state_min(), no_state_max());
    REQUIRE(eset.has_value());

    // Expected: alpha = min(1, 1) / (0.5^2 + 0.5^2) = 1 / 0.5 = 2.0
    CHECK_THAT(eset->alpha, WithinAbs(2.0, 1e-10));
}

TEST_CASE("inscribed-box vertices lie inside the terminal ellipsoid", "[terminal_set]")
{
    // Non-axis-aligned P so the box axes are the eigenvectors of P, not the coordinate axes.
    Eigen::Matrix2d P;
    P << 2.0, 0.5, 0.5, 1.0;
    const double alpha = 1.5;
    ctrlpp::ellipsoidal_set<double, NX> s{.P = P, .alpha = alpha};

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(P);
    auto V = eig.eigenvectors();
    auto D = eig.eigenvalues();

    // Inscribed-box axis bounds in eigen-coordinates: |y_i| <= sqrt(alpha / (NX * lambda_i)).
    Eigen::Vector2d bound;
    for(int i = 0; i < 2; ++i)
        bound(i) = std::sqrt(alpha / (static_cast<double>(NX) * D(i)));

    // Every box vertex x = V y (each y_i at +/- bound_i) must satisfy x^T P x <= alpha.
    for(int sx = -1; sx <= 1; sx += 2)
        for(int sy = -1; sy <= 1; sy += 2)
        {
            Eigen::Vector2d y{sx * bound(0), sy * bound(1)};
            Eigen::Vector2d x = V * y;
            double quad = x.transpose() * P * x;
            CHECK(quad <= alpha + 1e-10);
        }
    (void)s;
}

TEST_CASE("ellipsoidal terminal set encodes exactly NX rows", "[terminal_set]")
{
    ctrlpp::ellipsoidal_set<double, NX> s{.P = Eigen::Matrix2d::Identity(), .alpha = 1.0};
    std::optional<ctrlpp::terminal_set<double, NX>> tset{s};
    CHECK(ctrlpp::detail::terminal_constraint_rows<double, NX>(tset) == static_cast<int>(NX));
}

TEST_CASE("compute_ellipsoidal_set rejects a u-range that does not straddle zero", "[terminal_set]")
{
    Eigen::Matrix2d P = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 1, 2> K;
    K << 0.5, 0.5;
    Eigen::Matrix<double, 1, 1> u_min;
    u_min << 0.5; // u = 0 is not interior to [0.5, 2.0]
    Eigen::Matrix<double, 1, 1> u_max;
    u_max << 2.0;

    auto eset = ctrlpp::compute_ellipsoidal_set<double, NX, NU>(P, K, u_min, u_max, no_state_min(), no_state_max());
    REQUIRE_FALSE(eset.has_value());
    CHECK(eset.error() == ctrlpp::terminal_set_error::input_zero_not_interior);
}

TEST_CASE("state faces cap alpha below the input-face bound", "[terminal_set]")
{
    // Small gain and huge input limits give a large input-face alpha; the state box then binds.
    Eigen::Matrix2d P = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 1, 2> K;
    K << 0.1, 0.0;
    Eigen::Matrix<double, 1, 1> u_min;
    u_min << -10.0;
    Eigen::Matrix<double, 1, 1> u_max;
    u_max << 10.0;
    Eigen::Vector2d x_min{-1.0, -1.0};
    Eigen::Vector2d x_max{1.0, 1.0};

    auto eset = ctrlpp::compute_ellipsoidal_set<double, NX, NU>(P, K, u_min, u_max, x_min, x_max);
    REQUIRE(eset.has_value());

    // Input face: min(100, 100) / 0.01 = 10000. State face (P = I): min(1, 1) / 1 = 1.
    CHECK_THAT(eset->alpha, WithinAbs(1.0, 1e-10));
}

TEST_CASE("degenerate terminal set returns empty_terminal_set", "[terminal_set]")
{
    // Zero gain deactivates the input face; unbounded state box leaves alpha at +inf.
    Eigen::Matrix2d P = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 1, 2> K;
    K << 0.0, 0.0;
    Eigen::Matrix<double, 1, 1> u_min;
    u_min << -1.0;
    Eigen::Matrix<double, 1, 1> u_max;
    u_max << 1.0;

    auto eset = ctrlpp::compute_ellipsoidal_set<double, NX, NU>(P, K, u_min, u_max, no_state_min(), no_state_max());
    REQUIRE_FALSE(eset.has_value());
    CHECK(eset.error() == ctrlpp::terminal_set_error::empty_terminal_set);
}

TEST_CASE("filter_redundant_halfplanes flags the resource-cap truncation", "[terminal_set]")
{
    // Build more distinct, non-redundant halfplane directions than the resource cap allows.
    constexpr int n = 600;
    const double two_pi = 2.0 * std::acos(-1.0);
    Eigen::Matrix<double, Eigen::Dynamic, 2> H(n, 2);
    Eigen::VectorXd h(n);
    for(int i = 0; i < n; ++i)
    {
        double theta = two_pi * static_cast<double>(i) / static_cast<double>(n);
        H(i, 0) = std::cos(theta);
        H(i, 1) = std::sin(theta);
        h(i) = 1.0;
    }

    auto filtered = ctrlpp::detail::filter_redundant_halfplanes<double, NX>(H, h, 1e-6);
    CHECK(filtered.truncated);
    CHECK(filtered.H.rows() <= 500);
}

TEST_CASE("compute_polytopic_invariant_set signals non-convergence", "[terminal_set]")
{
    auto sys = make_double_integrator();

    Eigen::Matrix<double, 4, 2> H_state;
    H_state << 1.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, -1.0;
    Eigen::Vector4d h_state;
    h_state << 5.0, 5.0, 5.0, 5.0;
    ctrlpp::polytopic_set<double, NX> state_constr{.H = H_state, .h = h_state};

    Eigen::Matrix<double, 2, 1> H_input;
    H_input << 1.0, -1.0;
    Eigen::Vector2d h_input;
    h_input << 1.0, 1.0;
    ctrlpp::polytopic_set<double, NU> input_constr{.H = H_input, .h = h_input};

    // One iteration cannot reach the fixed point: expect a signaled not_converged, not a silent set.
    auto result = ctrlpp::compute_polytopic_invariant_set<double, NX, NU>(sys.A, sys.B, state_constr, input_constr, 1);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::terminal_set_error::not_converged);
}

TEST_CASE("compute_polytopic_invariant_set on 2D system", "[terminal_set]")
{
    // Strictly contractive dynamics with small input action: the state box below is already
    // robustly control-invariant, so backward reachability converges to the box itself. (The
    // marginally stable double integrator does NOT converge within the halfplane budget for this
    // box; that non-convergence is exercised separately and was previously masked by a silent set.)
    Eigen::Matrix2d A;
    A << 0.5, 0.0, 0.0, 0.5;
    Eigen::Vector2d B;
    B << 0.1, 0.1;

    // Box state constraints: |x_i| <= 5
    Eigen::Matrix<double, 4, 2> H_state;
    H_state << 1.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, -1.0;
    Eigen::Vector4d h_state;
    h_state << 5.0, 5.0, 5.0, 5.0;
    ctrlpp::polytopic_set<double, NX> state_constr{.H = H_state, .h = h_state};

    // Box input constraints: |u| <= 1
    Eigen::Matrix<double, 2, 1> H_input;
    H_input << 1.0, -1.0;
    Eigen::Vector2d h_input;
    h_input << 1.0, 1.0;
    ctrlpp::polytopic_set<double, NU> input_constr{.H = H_input, .h = h_input};

    auto result = ctrlpp::compute_polytopic_invariant_set<double, NX, NU>(A, B, state_constr, input_constr, 50);

    REQUIRE(result.has_value());

    // Invariant set should be non-empty (has constraints)
    CHECK(result->H.rows() > 0);
    CHECK(result->h.size() == result->H.rows());

    // Origin should be inside the invariant set (H*0 <= h => all h >= 0)
    for(int i = 0; i < static_cast<int>(result->h.size()); ++i)
        CHECK(result->h(i) >= -1e-10);
}

TEST_CASE("MPC with terminal_ingredients integration", "[terminal_set][mpc]")
{
    auto sys = make_double_integrator();
    constexpr int N = 20;

    Eigen::Matrix2d Q = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 0.1;
    Eigen::Matrix<double, 1, 1> u_min;
    u_min << -5.0;
    Eigen::Matrix<double, 1, 1> u_max;
    u_max << 5.0;

    auto ti = ctrlpp::terminal_ingredients<double, NX, NU>(sys.A, sys.B, Q, R, u_min, u_max, no_state_min(), no_state_max());
    REQUIRE(ti.has_value());

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = N,
        .Q = Q,
        .R = R,
        .Qf = ti->Qf,
        .u_min = u_min,
        .u_max = u_max,
        .terminal_constraint_set = ctrlpp::terminal_set<double, NX>{ti->set},
    };

    auto controller_result = OsqpMpc::create(sys, cfg);
    REQUIRE(controller_result.has_value());
    auto& controller = *controller_result;

    // Closed-loop regulation: verify stability (small initial state for feasibility)
    Eigen::Vector2d x{0.1, 0.05};
    for(int step = 0; step < 50; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = sys.A * x + sys.B * u.value().input;
    }

    CHECK(x.norm() < 0.05);
}

TEST_CASE("MPC with polytopic terminal set", "[terminal_set][mpc]")
{
    auto sys = make_double_integrator();
    constexpr int N = 20;

    // Terminal box: |x_i| <= 1.0
    Eigen::Matrix<double, 4, 2> H_term;
    H_term << 1.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, -1.0;
    Eigen::Vector4d h_term;
    h_term << 1.0, 1.0, 1.0, 1.0;

    ctrlpp::polytopic_set<double, NX> pset{.H = H_term, .h = h_term};

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = N,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        .u_min = (Eigen::Matrix<double, 1, 1>() << -5.0).finished(),
        .u_max = (Eigen::Matrix<double, 1, 1>() << 5.0).finished(),
        .terminal_constraint_set = ctrlpp::terminal_set<double, NX>{pset},
    };

    auto controller_result = OsqpMpc::create(sys, cfg);
    REQUIRE(controller_result.has_value());
    auto& controller = *controller_result;

    Eigen::Vector2d x{0.5, 0.1};
    auto u = controller.solve(x);
    auto diag = controller.diagnostics();
    INFO("Solve status: " << static_cast<int>(diag.status));
    REQUIRE(u.has_value());

    // Extract trajectory and check terminal state satisfies Hx <= h
    auto traj = controller.trajectory();
    REQUIRE(traj.has_value());
    auto& [states, inputs] = *traj;
    auto x_N = states.back();
    Eigen::VectorXd Hx = H_term * x_N;
    for(int i = 0; i < 4; ++i)
        CHECK(Hx(i) <= h_term(i) + 0.02);
}

TEST_CASE("MPC backward compatibility without terminal_constraint_set", "[terminal_set][mpc]")
{
    auto sys = make_double_integrator();
    constexpr int N = 10;

    // Config without terminal set -- should produce same results as before
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
    CHECK(diag.status == ctrlpp::solve_status::optimal);
}

TEST_CASE("Per-state soft penalty", "[terminal_set][mpc]")
{
    auto sys = make_double_integrator();
    constexpr int N = 10;

    Eigen::Vector2d per_state;
    per_state << 1e5, 1e2; // position penalized more than velocity

    ctrlpp::mpc_config<double, NX, NU> cfg_per_state{
        .horizon = N,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        .x_min = (Eigen::Vector2d() << -2.0, -2.0).finished(),
        .x_max = (Eigen::Vector2d() << 2.0, 2.0).finished(),
        .soft_state_penalty = per_state,
    };

    ctrlpp::mpc_config<double, NX, NU> cfg_uniform{
        .horizon = N,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        .x_min = (Eigen::Vector2d() << -2.0, -2.0).finished(),
        .x_max = (Eigen::Vector2d() << 2.0, 2.0).finished(),
    };

    auto ctrl_ps_result = OsqpMpc::create(sys, cfg_per_state);
    REQUIRE(ctrl_ps_result.has_value());
    auto& ctrl_ps = *ctrl_ps_result;
    auto ctrl_uni_result = OsqpMpc::create(sys, cfg_uniform);
    REQUIRE(ctrl_uni_result.has_value());
    auto& ctrl_uni = *ctrl_uni_result;

    // Start outside bounds to trigger slack usage
    Eigen::Vector2d x0{5.0, 0.0};
    auto u_ps = ctrl_ps.solve(x0);
    auto u_uni = ctrl_uni.solve(x0);

    REQUIRE(u_ps.has_value());
    REQUIRE(u_uni.has_value());

    // Per-state and uniform should produce different costs
    CHECK(ctrl_ps.diagnostics().cost != ctrl_uni.diagnostics().cost);
}
