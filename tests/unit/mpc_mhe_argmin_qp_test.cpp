#include "ctrlpp/mpc/argmin_qp_solver.h"

#include <catch2/catch_test_macros.hpp>

// Catch2 exits non-zero from a binary with no registered test, so this case stays outside the gate.
TEST_CASE("mpc and mhe on argmin_qp_solver gate", "[mpc][mhe][argmin][qp]")
{
#if defined(CTRLPP_HAS_ARGMIN_QP)
    SUCCEED("argmin QP header present: the mpc and mhe cases below run");
#else
    SUCCEED("argmin QP header absent: the mpc and mhe cases are compiled out");
#endif
}

#if defined(CTRLPP_HAS_ARGMIN_QP)

#include "hardening_helpers.h"

#include "ctrlpp/mhe.h"
#include "ctrlpp/mpc.h"

#include <Eigen/Dense>

#include <cmath>
#include <cstddef>
#include <type_traits>

static_assert(std::is_move_constructible_v<ctrlpp::argmin_qp_solver>, "mpc and mhe hold their solver by value");

namespace
{

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr std::size_t NY = 1;
constexpr std::size_t N = 5;
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

struct linear_dynamics
{
    auto operator()(const ctrlpp::Vector<double, NX>& x, const ctrlpp::Vector<double, NU>& u) const -> ctrlpp::Vector<double, NX>
    {
        ctrlpp::Vector<double, NX> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    }

    auto jacobian_x(const ctrlpp::Vector<double, NX>& /*x*/, const ctrlpp::Vector<double, NU>& /*u*/) const -> ctrlpp::Matrix<double, NX, NX>
    {
        ctrlpp::Matrix<double, NX, NX> F;
        F << 1.0, dt, 0.0, 1.0;
        return F;
    }

    auto jacobian_u(const ctrlpp::Vector<double, NX>& /*x*/, const ctrlpp::Vector<double, NU>& /*u*/) const -> ctrlpp::Matrix<double, NX, NU>
    {
        ctrlpp::Matrix<double, NX, NU> G;
        G << 0.5 * dt * dt, dt;
        return G;
    }
};

struct position_measurement
{
    auto operator()(const ctrlpp::Vector<double, NX>& x) const -> ctrlpp::Vector<double, NY>
    {
        return (ctrlpp::Vector<double, NY>() << x(0)).finished();
    }

    auto jacobian(const ctrlpp::Vector<double, NX>& /*x*/) const -> ctrlpp::Matrix<double, NY, NX>
    {
        return (ctrlpp::Matrix<double, NY, NX>() << 1.0, 0.0).finished();
    }
};

using argmin_mpc = ctrlpp::mpc<double, NX, NU, ctrlpp::argmin_qp_solver>;
using argmin_mhe = ctrlpp::mhe<double, NX, NU, NY, N, ctrlpp::argmin_qp_solver, linear_dynamics, position_measurement>;

}

TEST_CASE("mpc on argmin_qp_solver regulates the double integrator", "[mpc][argmin][qp]")
{
    const auto sys = make_double_integrator();
    const ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = 10,
        .Q = Eigen::Matrix2d::Identity(),
        .R = (Eigen::Matrix<double, 1, 1>() << 0.1).finished(),
    };

    auto controller = ctrlpp::test::constructed(argmin_mpc::create(sys, cfg, ctrlpp::argmin_qp_solver{ctrlpp::qp_preset::speed}));

    Eigen::Vector2d x{1.0, 0.0};
    for(std::size_t step = 0; step < 50; ++step)
    {
        const auto u = controller.solve(x);
        REQUIRE(u.has_value());
        x = sys.A * x + sys.B * u->input;
    }
    CHECK(x.norm() < 0.1);
}

TEST_CASE("mhe on argmin_qp_solver solves its window after warm-up", "[mhe][argmin][qp]")
{
    ctrlpp::mhe_config<double, NX, NU, NY, N> cfg;
    cfg.Q = ctrlpp::Matrix<double, NX, NX>::Identity() * 0.01;
    cfg.R = ctrlpp::Matrix<double, NY, NY>::Identity() * 0.1;
    cfg.P0 = ctrlpp::Matrix<double, NX, NX>::Identity() * 10.0;

    auto estimator = ctrlpp::test::constructed(argmin_mhe::create(linear_dynamics{}, position_measurement{}, cfg));

    ctrlpp::Vector<double, NX> x_true{0.0, 1.0};
    const ctrlpp::Vector<double, NU> u = ctrlpp::Vector<double, NU>::Zero();
    for(std::size_t i = 0; i < N + 10; ++i)
    {
        x_true(0) += dt * x_true(1);
        estimator.predict(u);
        REQUIRE(estimator.update((ctrlpp::Vector<double, NY>() << x_true(0)).finished()).has_value());

        // The window QP first runs on the update after the one that fills the N-sample window.
        const auto& diag = estimator.diagnostics();
        if(i < N)
        {
            CHECK(diag.used_ekf_fallback);
            continue;
        }
        CHECK_FALSE(diag.used_ekf_fallback);
        CHECK(diag.status == ctrlpp::solve_status::optimal);
    }

    const auto& est = estimator.state();
    CHECK(std::isfinite(est(0)));
    CHECK(std::isfinite(est(1)));
}

#endif
