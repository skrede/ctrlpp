// Construction-time horizon domain for the two runtime-horizon controllers.
//
// The prediction horizon is the one runtime quantity that scales every derived
// dimension of the posed problem, and it is caller-supplied. These cases pin the
// boundary of its accepted domain from both sides and assert the specific typed
// error on each rejection, never merely that construction "did not crash".
//
// The stub solvers keep this translation unit backend-free, so it builds and
// runs in the default (-fno-exceptions) tree where the real backends are absent.

#include "ctrlpp/mpc.h"
#include "ctrlpp/nmpc.h"
#include "ctrlpp/model/state_space.h"
#include "ctrlpp/mpc/nlp_formulation.h"

#include "stub_qp_solver.h"
#include "stub_nlp_solver.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>

#include <limits>
#include <memory>
#include <cstddef>

namespace
{

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr double dt = 0.1;

using linear_controller = ctrlpp::mpc<double, NX, NU, ctrlpp_test::stub_qp_solver<double>>;

auto double_integrator = [](const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>& u) -> Eigen::Vector2d
{ return Eigen::Vector2d{x(0) + dt * x(1), x(1) + dt * u(0)}; };

using nonlinear_controller = ctrlpp::nmpc_dynamic<double, NX, NU, ctrlpp_test::stub_nlp_solver<double>, decltype(double_integrator)>;

auto make_system() -> ctrlpp::discrete_state_space<double, NX, NU, NX>
{
    ctrlpp::discrete_state_space<double, NX, NU, NX> sys;
    sys.A << 1.0, dt, 0.0, 1.0;
    sys.B << 0.0, dt;
    sys.C = Eigen::Matrix2d::Identity();
    sys.D = Eigen::Matrix<double, 2, 1>::Zero();
    return sys;
}

auto linear_config(int horizon) -> ctrlpp::mpc_config<double, NX, NU>
{
    ctrlpp::mpc_config<double, NX, NU> cfg;
    cfg.horizon = horizon;
    return cfg;
}

auto nonlinear_config(int horizon) -> ctrlpp::nmpc_config<double, NX, NU>
{
    ctrlpp::nmpc_config<double, NX, NU> cfg;
    cfg.horizon = horizon;
    return cfg;
}

// The bounds below restate the derivations documented on the two factories, so a
// change to either derivation without a matching change here is caught rather
// than absorbed. Both are exact representability conditions on int, computed
// from the per-step dimension contribution; neither is a chosen ceiling.
constexpr int linear_horizon_bound = (std::numeric_limits<int>::max() - static_cast<int>(NX)) / (2 * static_cast<int>(NX) + 2 * static_cast<int>(NU));
constexpr int nonlinear_horizon_bound = (std::numeric_limits<int>::max() - static_cast<int>(NX)) / (static_cast<int>(NX) + 2 * static_cast<int>(NU));

}

TEST_CASE("Linear MPC rejects a negative horizon at construction", "[mpc][horizon][hardening]")
{
    auto result = linear_controller::try_create(make_system(), linear_config(-1), ctrlpp_test::stub_qp_solver<double>{});

    REQUIRE(!result.has_value());
    REQUIRE(result.error() == ctrlpp::controller_construction_error::non_positive_horizon);
}

TEST_CASE("Linear MPC rejects a zero horizon at construction", "[mpc][horizon][hardening]")
{
    auto result = linear_controller::try_create(make_system(), linear_config(0), ctrlpp_test::stub_qp_solver<double>{});

    REQUIRE(!result.has_value());
    REQUIRE(result.error() == ctrlpp::controller_construction_error::non_positive_horizon);
}

TEST_CASE("Linear MPC accepts a horizon of one and solves", "[mpc][horizon][hardening]")
{
    auto result = linear_controller::try_create(make_system(), linear_config(1), ctrlpp_test::stub_qp_solver<double>{});

    REQUIRE(result.has_value());

    auto solved = result->solve(Eigen::Vector2d{1.0, 0.0});
    REQUIRE(solved.has_value());
    REQUIRE(solved->status == ctrlpp::solve_result_status::converged);

    auto traj = result->trajectory();
    REQUIRE(traj.has_value());
    REQUIRE(traj->first.size() == 2);
    REQUIRE(traj->second.size() == 1);
}

TEST_CASE("Linear MPC rejects a horizon whose derived dimensions would overflow", "[mpc][horizon][hardening]")
{
    // The first horizon past the representability bound, and the largest value
    // the field can hold, are both rejected as overflow rather than as a
    // non-positive horizon. Constructing AT the bound is not attempted: it is a
    // valid horizon, so it would genuinely try to reserve a decision vector of
    // roughly two billion entries.
    auto just_past = linear_controller::try_create(make_system(), linear_config(linear_horizon_bound + 1), ctrlpp_test::stub_qp_solver<double>{});
    REQUIRE(!just_past.has_value());
    REQUIRE(just_past.error() == ctrlpp::controller_construction_error::horizon_overflow);

    auto largest = linear_controller::try_create(make_system(), linear_config(std::numeric_limits<int>::max()), ctrlpp_test::stub_qp_solver<double>{});
    REQUIRE(!largest.has_value());
    REQUIRE(largest.error() == ctrlpp::controller_construction_error::horizon_overflow);
}

TEST_CASE("Nonlinear MPC rejects a negative horizon at construction", "[nmpc][horizon][hardening]")
{
    auto result = nonlinear_controller::try_create(double_integrator, nonlinear_config(-1), ctrlpp_test::stub_nlp_solver<double>{});

    REQUIRE(!result.has_value());
    REQUIRE(result.error() == ctrlpp::controller_construction_error::non_positive_horizon);
}

TEST_CASE("Nonlinear MPC rejects a zero horizon at construction", "[nmpc][horizon][hardening]")
{
    auto result = nonlinear_controller::try_create(double_integrator, nonlinear_config(0), ctrlpp_test::stub_nlp_solver<double>{});

    REQUIRE(!result.has_value());
    REQUIRE(result.error() == ctrlpp::controller_construction_error::non_positive_horizon);
}

TEST_CASE("Nonlinear MPC accepts a horizon of one and solves", "[nmpc][horizon][hardening]")
{
    auto result = nonlinear_controller::try_create(double_integrator, nonlinear_config(1), ctrlpp_test::stub_nlp_solver<double>{});

    REQUIRE(result.has_value());

    auto solved = result->solve(Eigen::Vector2d{1.0, 0.0});
    REQUIRE(solved.has_value());
    REQUIRE(solved->status == ctrlpp::solve_result_status::converged);

    auto traj = result->trajectory();
    REQUIRE(traj.has_value());
    REQUIRE(traj->first.size() == 2);
    REQUIRE(traj->second.size() == 1);
}

TEST_CASE("Nonlinear MPC rejects a horizon whose derived dimensions would overflow", "[nmpc][horizon][hardening]")
{
    auto just_past = nonlinear_controller::try_create(double_integrator, nonlinear_config(nonlinear_horizon_bound + 1), ctrlpp_test::stub_nlp_solver<double>{});
    REQUIRE(!just_past.has_value());
    REQUIRE(just_past.error() == ctrlpp::controller_construction_error::horizon_overflow);

    auto largest = nonlinear_controller::try_create(double_integrator, nonlinear_config(std::numeric_limits<int>::max()), ctrlpp_test::stub_nlp_solver<double>{});
    REQUIRE(!largest.has_value());
    REQUIRE(largest.error() == ctrlpp::controller_construction_error::horizon_overflow);
}

TEST_CASE("Compile-time-horizon NLP formulation rejects a disagreeing runtime horizon", "[nmpc][horizon][hardening]")
{
    constexpr std::size_t NH = 4;

    auto state = std::make_shared<ctrlpp::nmpc_formulation_state<double, NX, NU>>();
    state->x_ref.assign(NH + 1, Eigen::Vector2d::Zero());

    SECTION("a runtime horizon below the compile-time one")
    {
        auto built = ctrlpp::detail::build_nmpc_problem_static<double, NX, NU, NH>(double_integrator, nonlinear_config(static_cast<int>(NH) - 1), state);
        REQUIRE(!built.has_value());
        REQUIRE(built.error() == ctrlpp::nlp_formulation_error::horizon_mismatch);
    }

    SECTION("a runtime horizon above the compile-time one")
    {
        auto built = ctrlpp::detail::build_nmpc_problem_static<double, NX, NU, NH>(double_integrator, nonlinear_config(static_cast<int>(NH) + 1), state);
        REQUIRE(!built.has_value());
        REQUIRE(built.error() == ctrlpp::nlp_formulation_error::horizon_mismatch);
    }

    SECTION("a non-positive runtime horizon")
    {
        auto built = ctrlpp::detail::build_nmpc_problem_static<double, NX, NU, NH>(double_integrator, nonlinear_config(0), state);
        REQUIRE(!built.has_value());
        REQUIRE(built.error() == ctrlpp::nlp_formulation_error::horizon_mismatch);
    }

    SECTION("the agreeing runtime horizon is accepted")
    {
        auto built = ctrlpp::detail::build_nmpc_problem_static<double, NX, NU, NH>(double_integrator, nonlinear_config(static_cast<int>(NH)), state);
        REQUIRE(built.has_value());
        REQUIRE(built->n_vars == static_cast<int>((NH + 1) * NX + NH * NU));
    }
}

TEST_CASE("Compile-time-horizon controller reports setup_incomplete on a disagreeing runtime horizon", "[nmpc][horizon][hardening]")
{
    constexpr std::size_t NH = 4;
    using static_controller = ctrlpp::nmpc_static<double, NX, NU, NH, ctrlpp_test::stub_nlp_solver<double>, decltype(double_integrator)>;

    static_controller controller{double_integrator, nonlinear_config(static_cast<int>(NH) + 1)};

    auto solved = controller.solve(Eigen::Vector2d{1.0, 0.0});
    REQUIRE(!solved.has_value());
    REQUIRE(solved.error() == ctrlpp::solver_error::setup_incomplete);
}
