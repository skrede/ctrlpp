#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/invariant_set.h"
#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/mpc/terminal_set.h"

#include "ctrlpp/dare.h"
#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>

namespace {

using Catch::Matchers::WithinAbs;

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr double dt = 0.1;

auto make_double_integrator() -> ctrlpp::discrete_state_space<double, NX, NU, NX>
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

using OsqpMpc = ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver>;

}

TEST_CASE("terminal_ingredients on stable double integrator", "[terminal_set]") {
    auto sys = make_double_integrator();
    Eigen::Matrix2d Q = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 0.1;
    Eigen::Matrix<double, 1, 1> u_min;
    u_min << -1.0;
    Eigen::Matrix<double, 1, 1> u_max;
    u_max << 1.0;

    auto result = ctrlpp::terminal_ingredients<double, NX, NU>(sys.A, sys.B, Q, R, u_min, u_max);

    REQUIRE(result.has_value());

    // Qf should be positive semi-definite
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(result->Qf);
    for (int i = 0; i < 2; ++i)
        CHECK(eig.eigenvalues()(i) >= -1e-10);

    // alpha should be positive and finite
    CHECK(result->set.alpha > 0.0);
    CHECK(std::isfinite(result->set.alpha));
}

TEST_CASE("compute_ellipsoidal_set with known DARE solution", "[terminal_set]") {
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

    auto eset = ctrlpp::compute_ellipsoidal_set<double, NX, NU>(P, K, u_min, u_max);

    // Expected: alpha = min(1, 1) / (0.5^2 + 0.5^2) = 1 / 0.5 = 2.0
    CHECK_THAT(eset.alpha, WithinAbs(2.0, 1e-10));
}

TEST_CASE("compute_polytopic_invariant_set on 2D system", "[terminal_set]") {
    auto sys = make_double_integrator();

    // Box state constraints: |x_i| <= 5
    Eigen::Matrix<double, 4, 2> H_state;
    H_state << 1.0, 0.0,
              -1.0, 0.0,
               0.0, 1.0,
               0.0,-1.0;
    Eigen::Vector4d h_state;
    h_state << 5.0, 5.0, 5.0, 5.0;
    ctrlpp::polytopic_set<double, NX> state_constr{.H = H_state, .h = h_state};

    // Box input constraints: |u| <= 1
    Eigen::Matrix<double, 2, 1> H_input;
    H_input << 1.0, -1.0;
    Eigen::Vector2d h_input;
    h_input << 1.0, 1.0;
    ctrlpp::polytopic_set<double, NU> input_constr{.H = H_input, .h = h_input};

    auto result = ctrlpp::compute_polytopic_invariant_set<double, NX, NU>(
        sys.A, sys.B, state_constr, input_constr, 50);

    REQUIRE(result.has_value());

    // Invariant set should be non-empty (has constraints)
    CHECK(result->H.rows() > 0);
    CHECK(result->h.size() == result->H.rows());

    // Origin should be inside the invariant set (H*0 <= h => all h >= 0)
    for (int i = 0; i < static_cast<int>(result->h.size()); ++i)
        CHECK(result->h(i) >= -1e-10);
}

TEST_CASE("MPC with terminal_ingredients integration", "[terminal_set][mpc]") {
    auto sys = make_double_integrator();
    constexpr int N = 20;

    Eigen::Matrix2d Q = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 0.1;
    Eigen::Matrix<double, 1, 1> u_min;
    u_min << -5.0;
    Eigen::Matrix<double, 1, 1> u_max;
    u_max << 5.0;

    auto ti = ctrlpp::terminal_ingredients<double, NX, NU>(sys.A, sys.B, Q, R, u_min, u_max);
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

    OsqpMpc controller(sys, cfg);

    // Closed-loop regulation: verify stability (small initial state for feasibility)
    Eigen::Vector2d x{0.1, 0.05};
    for (int step = 0; step < 50; ++step) {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = sys.A * x + sys.B * u.value();
    }

    CHECK(x.norm() < 0.05);
}

TEST_CASE("MPC with polytopic terminal set", "[terminal_set][mpc]") {
    auto sys = make_double_integrator();
    constexpr int N = 20;

    // Terminal box: |x_i| <= 1.0
    Eigen::Matrix<double, 4, 2> H_term;
    H_term << 1.0, 0.0,
             -1.0, 0.0,
              0.0, 1.0,
              0.0,-1.0;
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

    OsqpMpc controller(sys, cfg);

    Eigen::Vector2d x{0.5, 0.1};
    auto u = controller.solve(x);
    auto diag = controller.diagnostics();
    INFO("Solve status: " << static_cast<int>(diag.status));
    REQUIRE(u.has_value());

    // Extract trajectory and check terminal state satisfies Hx <= h
    auto [states, inputs] = controller.trajectory();
    auto x_N = states.back();
    Eigen::VectorXd Hx = H_term * x_N;
    for (int i = 0; i < 4; ++i)
        CHECK(Hx(i) <= h_term(i) + 0.02);
}

TEST_CASE("MPC backward compatibility without terminal_constraint_set", "[terminal_set][mpc]") {
    auto sys = make_double_integrator();
    constexpr int N = 10;

    // Config without terminal set -- should produce same results as before
    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = N,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
    };

    OsqpMpc controller(sys, cfg);

    Eigen::Vector2d x0{1.0, 0.0};
    auto result = controller.solve(x0);
    REQUIRE(result.has_value());

    auto diag = controller.diagnostics();
    CHECK(diag.status == ctrlpp::solve_status::optimal);
}

TEST_CASE("Per-state soft penalty", "[terminal_set][mpc]") {
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

    OsqpMpc ctrl_ps(sys, cfg_per_state);
    OsqpMpc ctrl_uni(sys, cfg_uniform);

    // Start outside bounds to trigger slack usage
    Eigen::Vector2d x0{5.0, 0.0};
    auto u_ps = ctrl_ps.solve(x0);
    auto u_uni = ctrl_uni.solve(x0);

    REQUIRE(u_ps.has_value());
    REQUIRE(u_uni.has_value());

    // Per-state and uniform should produce different costs
    CHECK(ctrl_ps.diagnostics().cost != ctrl_uni.diagnostics().cost);
}
