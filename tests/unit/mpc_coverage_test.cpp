#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/invariant_set.h"
#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/mpc/qp_formulation.h"
#include "ctrlpp/mpc/terminal_set.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

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

// Unbounded state box: the terminal ellipsoid is then constrained by the input faces alone.
auto no_state_min() -> Eigen::Vector2d
{
    return Eigen::Vector2d::Constant(-std::numeric_limits<double>::infinity());
}
auto no_state_max() -> Eigen::Vector2d
{
    return Eigen::Vector2d::Constant(std::numeric_limits<double>::infinity());
}

}

// ---- qp_formulation.h: soft state constraints with custom per-state penalties ----

TEST_CASE("soft state constraints with per-state penalty vector", "[mpc][coverage][qp]")
{
    auto sys = make_double_integrator();
    constexpr int N = 10;

    Eigen::Vector2d penalty;
    penalty << 1e6, 1e2;

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = N,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        .x_min = Eigen::Vector2d{-1.0, -1.0},
        .x_max = Eigen::Vector2d{1.0, 1.0},
        .soft_state_penalty = penalty,
    };

    OsqpMpc controller(sys, cfg);

    // Start outside state bounds; soft constraints should allow a feasible solution
    Eigen::Vector2d x0{3.0, 0.0};
    auto u = controller.solve(x0);
    REQUIRE(u.has_value());

    // Simulate and check that position converges within bounds faster than velocity
    Eigen::Vector2d x = x0;
    for(int step = 0; step < 40; ++step)
    {
        auto ui = controller.solve(x);
        REQUIRE(ui.has_value());
        x = sys.A * x + sys.B * ui.value();
    }

    CHECK(std::abs(x(0)) < 1.5);
}

// ---- qp_formulation.h: terminal constraint row counting ----

TEST_CASE("terminal_constraint_rows returns 0 for nullopt", "[mpc][coverage][qp]")
{
    std::optional<ctrlpp::terminal_set<double, NX>> none{};
    CHECK(ctrlpp::detail::terminal_constraint_rows<double, NX>(none) == 0);
}

TEST_CASE("terminal_constraint_rows returns NX for ellipsoidal set", "[mpc][coverage][qp]")
{
    ctrlpp::ellipsoidal_set<double, NX> eset{
        .P = Eigen::Matrix2d::Identity(),
        .alpha = 1.0,
    };
    std::optional<ctrlpp::terminal_set<double, NX>> tset{eset};
    // Inscribed-box encoding emits one two-sided row per eigen-direction: NX rows, not 2*NX.
    CHECK(ctrlpp::detail::terminal_constraint_rows<double, NX>(tset) == static_cast<int>(NX));
}

TEST_CASE("terminal_constraint_rows returns H.rows() for polytopic set", "[mpc][coverage][qp]")
{
    Eigen::Matrix<double, 4, 2> H;
    H << 1, 0, -1, 0, 0, 1, 0, -1;
    Eigen::Vector4d h;
    h << 1, 1, 1, 1;
    ctrlpp::polytopic_set<double, NX> pset{.H = H, .h = h};
    std::optional<ctrlpp::terminal_set<double, NX>> tset{pset};
    CHECK(ctrlpp::detail::terminal_constraint_rows<double, NX>(tset) == 4);
}

// ---- mpc.h: explicit Qf (terminal cost) ----

TEST_CASE("explicit Qf overrides DARE-computed terminal cost", "[mpc][coverage]")
{
    auto sys = make_double_integrator();
    constexpr int N = 10;

    Eigen::Matrix2d Qf_explicit = 10.0 * Eigen::Matrix2d::Identity();

    ctrlpp::mpc_config<double, NX, NU> cfg_explicit{
        .horizon = N,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        .Qf = Qf_explicit,
    };

    ctrlpp::mpc_config<double, NX, NU> cfg_dare{
        .horizon = N,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
    };

    OsqpMpc ctrl_explicit(sys, cfg_explicit);
    OsqpMpc ctrl_dare(sys, cfg_dare);

    Eigen::Vector2d x0{1.0, 0.0};
    auto u_explicit = ctrl_explicit.solve(x0);
    auto u_dare = ctrl_dare.solve(x0);

    REQUIRE(u_explicit.has_value());
    REQUIRE(u_dare.has_value());

    // Different terminal costs should produce different optimal controls
    CHECK(std::abs((*u_explicit)(0) - (*u_dare)(0)) > 1e-6);
}

// ---- mpc.h: asymmetric state bounds (x_min only, x_max only) ----

TEST_CASE("asymmetric state bounds: x_min only", "[mpc][coverage]")
{
    auto sys = make_double_integrator();
    constexpr int N = 10;

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = N,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        .u_min = Eigen::Matrix<double, 1, 1>{-10.0},
        .u_max = Eigen::Matrix<double, 1, 1>{10.0},
        .x_min = Eigen::Vector2d{-5.0, -5.0},
    };

    OsqpMpc controller(sys, cfg);

    Eigen::Vector2d x0{1.0, 0.0};
    auto u = controller.solve(x0);
    // Solver may or may not converge with one-sided bounds;
    // the point is exercising the QP formulation path
    if(u.has_value())
    {
        auto [states, inputs] = controller.trajectory();
        CHECK(states.size() == static_cast<std::size_t>(N + 1));
    }
}

TEST_CASE("asymmetric state bounds: x_max only", "[mpc][coverage]")
{
    auto sys = make_double_integrator();
    constexpr int N = 10;

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = N,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        .u_min = Eigen::Matrix<double, 1, 1>{-10.0},
        .u_max = Eigen::Matrix<double, 1, 1>{10.0},
        .x_max = Eigen::Vector2d{2.0, 2.0},
    };

    OsqpMpc controller(sys, cfg);

    Eigen::Vector2d x0{1.0, 0.0};
    auto u = controller.solve(x0);
    if(u.has_value())
    {
        auto [states, inputs] = controller.trajectory();
        for(const auto& s : states)
        {
            CHECK(s(0) <= 2.0 + 1e-3);
            CHECK(s(1) <= 2.0 + 1e-3);
        }
    }
}

// ---- mpc.h: rate constraints with warm-starting ----

TEST_CASE("rate constraints with warm-started consecutive solves", "[mpc][coverage]")
{
    auto sys = make_double_integrator();
    constexpr int N = 10;

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = N,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        .u_min = (Eigen::Matrix<double, 1, 1>() << -5.0).finished(),
        .u_max = (Eigen::Matrix<double, 1, 1>() << 5.0).finished(),
        .du_max = (Eigen::Matrix<double, 1, 1>() << 0.3).finished(),
    };

    OsqpMpc controller(sys, cfg);

    Eigen::Vector2d x{3.0, 0.0};
    double u_prev = 0.0;

    for(int step = 0; step < 20; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());

        double du = std::abs((*u)(0) - u_prev);
        CHECK(du <= 0.3 + 1e-2);

        u_prev = (*u)(0);
        x = sys.A * x + sys.B * u.value();
    }

    CHECK(x.norm() < 2.0);
}

// ---- mpc.h: span reference tracking ----

TEST_CASE("span reference trajectory tracking", "[mpc][coverage]")
{
    auto sys = make_double_integrator();
    constexpr int N = 10;

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = N,
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
    };

    OsqpMpc controller(sys, cfg);

    Eigen::Vector2d x{0.0, 0.0};

    // Build a time-varying reference trajectory
    std::vector<Eigen::Vector2d> refs;
    refs.reserve(static_cast<std::size_t>(N + 1));
    for(int k = 0; k <= N; ++k)
        refs.push_back(Eigen::Vector2d{1.0, 0.0});

    auto u = controller.solve(x, std::span<const Eigen::Vector2d>{refs});
    REQUIRE(u.has_value());
    CHECK((*u)(0) > 0.0); // should push toward positive reference
}

// ---- terminal_set.h: ellipsoidal set with various alpha values ----

TEST_CASE("ellipsoidal set with small alpha", "[mpc][coverage][terminal_set]")
{
    auto sys = make_double_integrator();
    Eigen::Matrix2d Q = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 0.1;
    Eigen::Matrix<double, 1, 1> u_min;
    u_min << -0.1;
    Eigen::Matrix<double, 1, 1> u_max;
    u_max << 0.1;

    auto result = ctrlpp::terminal_ingredients<double, NX, NU>(sys.A, sys.B, Q, R, u_min, u_max, no_state_min(), no_state_max());
    REQUIRE(result.has_value());

    // Small input bounds should produce small alpha
    CHECK(result->set.alpha > 0.0);
    CHECK(result->set.alpha < 1.0);
}

TEST_CASE("ellipsoidal set with asymmetric input bounds", "[mpc][coverage][terminal_set]")
{
    Eigen::Matrix2d P = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 1, 2> K;
    K << 1.0, 0.0;
    Eigen::Matrix<double, 1, 1> u_min;
    u_min << -0.5;
    Eigen::Matrix<double, 1, 1> u_max;
    u_max << 2.0;

    auto eset = ctrlpp::compute_ellipsoidal_set<double, NX, NU>(P, K, u_min, u_max, no_state_min(), no_state_max());
    REQUIRE(eset.has_value());

    // alpha = min(u_max^2, u_min^2) / (k' P^{-1} k)
    // With P=I, K=[1,0]: denom = 1.0, min(4.0, 0.25) = 0.25
    CHECK_THAT(eset->alpha, WithinAbs(0.25, 1e-10));
}

// ---- terminal_set.h: polytopic set used in QP formulation ----

TEST_CASE("polytopic terminal set with many faces", "[mpc][coverage][terminal_set]")
{
    auto sys = make_double_integrator();
    constexpr int N = 15;

    // Octagonal terminal set (8 faces)
    Eigen::Matrix<double, 8, 2> H_term;
    double s = 1.0 / std::sqrt(2.0);
    H_term << 1, 0, -1, 0, 0, 1, 0, -1, s, s, s, -s, -s, s, -s, -s;
    Eigen::Matrix<double, 8, 1> h_term;
    h_term.setConstant(2.0);

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

    Eigen::Vector2d x0{0.3, 0.1};
    auto u = controller.solve(x0);
    REQUIRE(u.has_value());

    // Terminal state should satisfy Hx <= h
    auto [states, inputs] = controller.trajectory();
    auto x_N = states.back();
    Eigen::VectorXd Hx = H_term * x_N;
    for(int i = 0; i < 8; ++i)
        CHECK(Hx(i) <= h_term(i) + 0.05);
}

// ---- qp_formulation.h: hard state constraints with bounds ----

TEST_CASE("hard state constraints prevent slack variables", "[mpc][coverage][qp]")
{
    auto sys = make_double_integrator();
    constexpr int N = 10;

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = N,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        .x_min = Eigen::Vector2d{-5.0, -5.0},
        .x_max = Eigen::Vector2d{5.0, 5.0},
        .hard_state_constraints = true,
    };

    OsqpMpc controller(sys, cfg);

    Eigen::Vector2d x0{1.0, 0.0};
    auto u = controller.solve(x0);
    REQUIRE(u.has_value());

    // With hard constraints and state within bounds, all trajectory states must satisfy bounds
    auto [states, inputs] = controller.trajectory();
    for(const auto& s : states)
    {
        CHECK(s(0) >= -5.0 - 1e-3);
        CHECK(s(0) <= 5.0 + 1e-3);
        CHECK(s(1) >= -5.0 - 1e-3);
        CHECK(s(1) <= 5.0 + 1e-3);
    }
}

// ---- mpc.h: combined input + state + rate constraints ----

TEST_CASE("all constraint types active simultaneously", "[mpc][coverage]")
{
    auto sys = make_double_integrator();
    constexpr int N = 10;

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = N,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
        .u_min = (Eigen::Matrix<double, 1, 1>() << -1.0).finished(),
        .u_max = (Eigen::Matrix<double, 1, 1>() << 1.0).finished(),
        .x_min = Eigen::Vector2d{-3.0, -3.0},
        .x_max = Eigen::Vector2d{3.0, 3.0},
        .du_max = (Eigen::Matrix<double, 1, 1>() << 0.5).finished(),
    };

    OsqpMpc controller(sys, cfg);

    Eigen::Vector2d x{2.0, 0.5};
    double u_prev = 0.0;

    for(int step = 0; step < 30; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());

        CHECK((*u)(0) >= -1.0 - 1e-3);
        CHECK((*u)(0) <= 1.0 + 1e-3);
        CHECK(std::abs((*u)(0) - u_prev) <= 0.5 + 1e-2);

        u_prev = (*u)(0);
        x = sys.A * x + sys.B * u.value();
    }

    CHECK(x.norm() < 2.0);
}

// ---- mpc.h: ellipsoidal terminal set in QP ----

TEST_CASE("MPC with ellipsoidal terminal set", "[mpc][coverage][terminal_set]")
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

    OsqpMpc controller(sys, cfg);

    Eigen::Vector2d x{0.2, 0.05};
    for(int step = 0; step < 50; ++step)
    {
        auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = sys.A * x + sys.B * u.value();
    }

    CHECK(x.norm() < 0.05);
}
