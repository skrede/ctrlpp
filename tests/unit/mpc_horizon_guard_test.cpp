// Two guards on the same data path: the horizon that sizes the posed problem,
// and the shape of the result the backend hands back.
//
// Construction-time horizon domain. The prediction horizon is the one runtime
// quantity that scales every derived dimension of the posed problem, and it is
// caller-supplied. These cases pin the boundary of its accepted domain from both
// sides and assert the specific typed error on each rejection, never merely that
// construction "did not crash".
//
// Result-shape domain. A valid horizon is not enough: the solve extracts fixed
// width slices out of whatever vector the backend returns, at offsets derived
// from the horizon, and today it does so on the strength of the reported status
// alone. A backend that reports an optimal status and returns a primal shorter
// than the decision dimension therefore reads past the end of that vector. The
// cases below drive exactly that, one per extraction site, and assert the
// specific observable each consumer reports.
//
// The stub solvers keep this translation unit backend-free, so it builds and
// runs in the default (-fno-exceptions) tree where the real backends are absent.

#include "ctrlpp/mhe.h"
#include "ctrlpp/mpc.h"
#include "ctrlpp/nmhe.h"
#include "ctrlpp/nmpc.h"
#include "ctrlpp/model/state_space.h"
#include "ctrlpp/mpc/nlp_formulation.h"

#include "hardening_helpers.h"
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

// The horizon the result-shape cases pose the problem at. It is a perfectly
// valid horizon: these cases are about the backend's answer, not the problem.
constexpr int shape_horizon = 5;

// Compile-time-horizon controller under a short-returning backend. The static
// path leaves its primal in a pre-sized buffer when the solver offers a
// write-into entry point and copies a by-value primal of unknown length when it
// does not; the stub offers only the latter, which is the copying shape.
template <ctrlpp_test::report_lengths Reported>
using compile_time_controller = ctrlpp::nmpc_static<double, NX, NU, 5, ctrlpp_test::stub_nlp_solver<double, Reported>, decltype(double_integrator)>;

// Estimator fixtures. Both estimators default-construct their solver member, so
// the reported result length is chosen as a template argument on the stub type.
constexpr std::size_t NY = 1;
constexpr std::size_t window = 3;

struct window_dynamics
{
    auto operator()(const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>& u) const -> Eigen::Vector2d { return Eigen::Vector2d{x(0) + dt * x(1), x(1) + dt * u(0)}; }
};

struct window_measurement
{
    auto operator()(const Eigen::Vector2d& x) const -> Eigen::Matrix<double, 1, 1> { return x.head<1>(); }
};

template <ctrlpp_test::report_lengths Reported>
using linear_estimator = ctrlpp::mhe<double, NX, NU, NY, window, ctrlpp_test::stub_qp_solver<double, Reported>, window_dynamics, window_measurement>;

template <ctrlpp_test::report_lengths Reported>
using nonlinear_estimator = ctrlpp::nmhe<double, NX, NU, NY, window, ctrlpp_test::stub_nlp_solver<double, Reported>, window_dynamics, window_measurement>;

// Fill the window and then take one more step, so the estimator leaves its
// warm-up branch (which reports the embedded filter's estimate directly) and
// actually runs a solve whose result is extracted.
template <typename Estimator>
void drive_past_warmup(Estimator& estimator)
{
    for(std::size_t k = 0; k <= window; ++k)
    {
        estimator.predict(Eigen::Matrix<double, 1, 1>::Zero());
        REQUIRE(estimator.update(Eigen::Matrix<double, 1, 1>::Zero()).has_value());
    }
}

}

TEST_CASE("Linear MPC rejects a negative horizon at construction", "[mpc][horizon][hardening]")
{
    auto result = linear_controller::create(make_system(), linear_config(-1), ctrlpp_test::stub_qp_solver<double>{});

    REQUIRE(!result.has_value());
    REQUIRE(result.error() == ctrlpp::controller_construction_error::non_positive_horizon);
}

TEST_CASE("Linear MPC rejects a zero horizon at construction", "[mpc][horizon][hardening]")
{
    auto result = linear_controller::create(make_system(), linear_config(0), ctrlpp_test::stub_qp_solver<double>{});

    REQUIRE(!result.has_value());
    REQUIRE(result.error() == ctrlpp::controller_construction_error::non_positive_horizon);
}

TEST_CASE("Linear MPC accepts a horizon of one and solves", "[mpc][horizon][hardening]")
{
    auto result = linear_controller::create(make_system(), linear_config(1), ctrlpp_test::stub_qp_solver<double>{});

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
    auto just_past = linear_controller::create(make_system(), linear_config(linear_horizon_bound + 1), ctrlpp_test::stub_qp_solver<double>{});
    REQUIRE(!just_past.has_value());
    REQUIRE(just_past.error() == ctrlpp::controller_construction_error::horizon_overflow);

    auto largest = linear_controller::create(make_system(), linear_config(std::numeric_limits<int>::max()), ctrlpp_test::stub_qp_solver<double>{});
    REQUIRE(!largest.has_value());
    REQUIRE(largest.error() == ctrlpp::controller_construction_error::horizon_overflow);
}

TEST_CASE("Nonlinear MPC rejects a negative horizon at construction", "[nmpc][horizon][hardening]")
{
    auto result = nonlinear_controller::create(double_integrator, nonlinear_config(-1), ctrlpp_test::stub_nlp_solver<double>{});

    REQUIRE(!result.has_value());
    REQUIRE(result.error() == ctrlpp::controller_construction_error::non_positive_horizon);
}

TEST_CASE("Nonlinear MPC rejects a zero horizon at construction", "[nmpc][horizon][hardening]")
{
    auto result = nonlinear_controller::create(double_integrator, nonlinear_config(0), ctrlpp_test::stub_nlp_solver<double>{});

    REQUIRE(!result.has_value());
    REQUIRE(result.error() == ctrlpp::controller_construction_error::non_positive_horizon);
}

TEST_CASE("Nonlinear MPC accepts a horizon of one and solves", "[nmpc][horizon][hardening]")
{
    auto result = nonlinear_controller::create(double_integrator, nonlinear_config(1), ctrlpp_test::stub_nlp_solver<double>{});

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
    auto just_past = nonlinear_controller::create(double_integrator, nonlinear_config(nonlinear_horizon_bound + 1), ctrlpp_test::stub_nlp_solver<double>{});
    REQUIRE(!just_past.has_value());
    REQUIRE(just_past.error() == ctrlpp::controller_construction_error::horizon_overflow);

    auto largest = nonlinear_controller::create(double_integrator, nonlinear_config(std::numeric_limits<int>::max()), ctrlpp_test::stub_nlp_solver<double>{});
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

// ---------------------------------------------------------------------------
// Result shape. The horizon is valid throughout this group; what varies is the
// length of the vector the backend returns while reporting an optimal status.
// ---------------------------------------------------------------------------

TEST_CASE("Linear MPC rejects a backend primal too short for the decision vector", "[mpc][result-shape][hardening]")
{
    SECTION("a primal one entry short of the decision dimension")
    {
        auto controller = linear_controller::create(make_system(), linear_config(shape_horizon), ctrlpp_test::stub_qp_solver<double>{ctrlpp_test::report_lengths::short_primal});
        REQUIRE(controller.has_value());

        auto solved = controller->solve(Eigen::Vector2d{1.0, 0.0});
        REQUIRE(!solved.has_value());
        REQUIRE(solved.error() == ctrlpp::solver_error::invalid_backend_result);
    }

    SECTION("a primal that is empty despite the reported optimal status")
    {
        // The degenerate end of the same defect, and the shape that first
        // exposed it: a controller with a valid multi-step horizon, handed a
        // primal far too short to hold the input block the extraction slices
        // out of it.
        auto controller = linear_controller::create(make_system(), linear_config(shape_horizon), ctrlpp_test::stub_qp_solver<double>{ctrlpp_test::report_lengths::empty});
        REQUIRE(controller.has_value());

        auto solved = controller->solve(Eigen::Vector2d{1.0, 0.0});
        REQUIRE(!solved.has_value());
        REQUIRE(solved.error() == ctrlpp::solver_error::invalid_backend_result);
    }
}

TEST_CASE("Linear MPC rejects a backend dual too short for the constraint rows", "[mpc][result-shape][hardening]")
{
    // The dual is not read by the extraction, it is stored and handed back to
    // the backend as the next warm start, so an undersized one is a deferred
    // overrun rather than an immediate one. It is rejected on the same branch.
    auto controller = linear_controller::create(make_system(), linear_config(shape_horizon), ctrlpp_test::stub_qp_solver<double>{ctrlpp_test::report_lengths::short_dual});
    REQUIRE(controller.has_value());

    auto solved = controller->solve(Eigen::Vector2d{1.0, 0.0});
    REQUIRE(!solved.has_value());
    REQUIRE(solved.error() == ctrlpp::solver_error::invalid_backend_result);
}

TEST_CASE("Linear MPC accepts a conforming backend result unchanged", "[mpc][result-shape][hardening]")
{
    auto controller = linear_controller::create(make_system(), linear_config(shape_horizon), ctrlpp_test::stub_qp_solver<double>{});
    REQUIRE(controller.has_value());

    auto solved = controller->solve(Eigen::Vector2d{1.0, 0.0});
    REQUIRE(solved.has_value());
    REQUIRE(solved->status == ctrlpp::solve_result_status::converged);

    auto traj = controller->trajectory();
    REQUIRE(traj.has_value());
    REQUIRE(traj->first.size() == static_cast<std::size_t>(shape_horizon) + 1);
    REQUIRE(traj->second.size() == static_cast<std::size_t>(shape_horizon));
}

TEST_CASE("Runtime-horizon nonlinear MPC rejects a backend primal too short for the problem", "[nmpc][result-shape][hardening]")
{
    SECTION("a primal one entry short of the problem dimension")
    {
        auto controller = nonlinear_controller::create(double_integrator, nonlinear_config(shape_horizon), ctrlpp_test::stub_nlp_solver<double>{ctrlpp_test::report_lengths::short_primal});
        REQUIRE(controller.has_value());

        auto solved = controller->solve(Eigen::Vector2d{1.0, 0.0});
        REQUIRE(!solved.has_value());
        REQUIRE(solved.error() == ctrlpp::solver_error::invalid_backend_result);
    }

    SECTION("a primal that is empty despite the reported optimal status")
    {
        auto controller = nonlinear_controller::create(double_integrator, nonlinear_config(shape_horizon), ctrlpp_test::stub_nlp_solver<double>{ctrlpp_test::report_lengths::empty});
        REQUIRE(controller.has_value());

        auto solved = controller->solve(Eigen::Vector2d{1.0, 0.0});
        REQUIRE(!solved.has_value());
        REQUIRE(solved.error() == ctrlpp::solver_error::invalid_backend_result);
    }
}

TEST_CASE("Runtime-horizon nonlinear MPC accepts a conforming backend result unchanged", "[nmpc][result-shape][hardening]")
{
    auto controller = nonlinear_controller::create(double_integrator, nonlinear_config(shape_horizon), ctrlpp_test::stub_nlp_solver<double>{});
    REQUIRE(controller.has_value());

    auto solved = controller->solve(Eigen::Vector2d{1.0, 0.0});
    REQUIRE(solved.has_value());
    REQUIRE(solved->status == ctrlpp::solve_result_status::converged);

    auto traj = controller->trajectory();
    REQUIRE(traj.has_value());
    REQUIRE(traj->first.size() == static_cast<std::size_t>(shape_horizon) + 1);
    REQUIRE(traj->second.size() == static_cast<std::size_t>(shape_horizon));
}

TEST_CASE("Compile-time-horizon nonlinear MPC rejects a backend primal too short for the problem", "[nmpc][result-shape][hardening]")
{
    SECTION("a primal one entry short of the compile-time dimension")
    {
        compile_time_controller<ctrlpp_test::report_lengths::short_primal> controller{double_integrator, nonlinear_config(shape_horizon)};

        auto solved = controller.solve(Eigen::Vector2d{1.0, 0.0});
        REQUIRE(!solved.has_value());
        REQUIRE(solved.error() == ctrlpp::solver_error::invalid_backend_result);
    }

    SECTION("a primal that is empty despite the reported optimal status")
    {
        compile_time_controller<ctrlpp_test::report_lengths::empty> controller{double_integrator, nonlinear_config(shape_horizon)};

        auto solved = controller.solve(Eigen::Vector2d{1.0, 0.0});
        REQUIRE(!solved.has_value());
        REQUIRE(solved.error() == ctrlpp::solver_error::invalid_backend_result);
    }
}

TEST_CASE("Compile-time-horizon nonlinear MPC accepts a conforming backend result unchanged", "[nmpc][result-shape][hardening]")
{
    compile_time_controller<ctrlpp_test::report_lengths::conforming> controller{double_integrator, nonlinear_config(shape_horizon)};

    auto solved = controller.solve(Eigen::Vector2d{1.0, 0.0});
    REQUIRE(solved.has_value());
    REQUIRE(solved->status == ctrlpp::solve_result_status::converged);
    REQUIRE(controller.last_solution().size() == decltype(controller)::problem_dimension);
}

TEST_CASE("Linear moving-horizon estimator falls back when the backend primal is too short", "[mhe][result-shape][hardening]")
{
    // The estimator's update returns nothing, so its failure channel is the
    // embedded-filter fallback plus the reported diagnostic status. Both are
    // asserted here; neither case asserts merely that the update returned.
    SECTION("a primal one entry short of the decision dimension")
    {
        auto estimator = ctrlpp::test::constructed(linear_estimator<ctrlpp_test::report_lengths::short_primal>::create(window_dynamics{}, window_measurement{}, ctrlpp::mhe_config<double, NX, NU, NY, window>{}));
        drive_past_warmup(estimator);

        REQUIRE(estimator.diagnostics().used_ekf_fallback);
        REQUIRE(estimator.diagnostics().status == ctrlpp::solve_status::invalid_backend_result);
    }

    SECTION("a primal that is empty despite the reported optimal status")
    {
        auto estimator = ctrlpp::test::constructed(linear_estimator<ctrlpp_test::report_lengths::empty>::create(window_dynamics{}, window_measurement{}, ctrlpp::mhe_config<double, NX, NU, NY, window>{}));
        drive_past_warmup(estimator);

        REQUIRE(estimator.diagnostics().used_ekf_fallback);
        REQUIRE(estimator.diagnostics().status == ctrlpp::solve_status::invalid_backend_result);
    }
}

TEST_CASE("Linear moving-horizon estimator falls back when the backend dual is too short", "[mhe][result-shape][hardening]")
{
    // State bounds are what give this problem constraint rows at all; without
    // them the constraint count is zero and no dual can be short of it.
    ctrlpp::mhe_config<double, NX, NU, NY, window> config;
    config.x_min = Eigen::Vector2d{-10.0, -10.0};
    config.x_max = Eigen::Vector2d{10.0, 10.0};

    auto estimator = ctrlpp::test::constructed(linear_estimator<ctrlpp_test::report_lengths::short_dual>::create(window_dynamics{}, window_measurement{}, config));
    drive_past_warmup(estimator);

    REQUIRE(estimator.diagnostics().used_ekf_fallback);
    REQUIRE(estimator.diagnostics().status == ctrlpp::solve_status::invalid_backend_result);
}

TEST_CASE("Linear moving-horizon estimator accepts a conforming backend result unchanged", "[mhe][result-shape][hardening]")
{
    auto estimator = ctrlpp::test::constructed(linear_estimator<ctrlpp_test::report_lengths::conforming>::create(window_dynamics{}, window_measurement{}, ctrlpp::mhe_config<double, NX, NU, NY, window>{}));
    drive_past_warmup(estimator);

    REQUIRE(!estimator.diagnostics().used_ekf_fallback);
    REQUIRE(estimator.diagnostics().status == ctrlpp::solve_status::optimal);
}

TEST_CASE("Nonlinear moving-horizon estimator falls back when the backend primal is too short", "[nmhe][result-shape][hardening]")
{
    SECTION("a primal one entry short of the decision dimension")
    {
        auto estimator = ctrlpp::test::constructed(nonlinear_estimator<ctrlpp_test::report_lengths::short_primal>::create(window_dynamics{}, window_measurement{}, ctrlpp::nmhe_config<double, NX, NU, NY, window>{}));
        drive_past_warmup(estimator);

        REQUIRE(estimator.diagnostics().used_ekf_fallback);
        REQUIRE(estimator.diagnostics().status == ctrlpp::solve_status::invalid_backend_result);
    }

    SECTION("a primal that is empty despite the reported optimal status")
    {
        auto estimator = ctrlpp::test::constructed(nonlinear_estimator<ctrlpp_test::report_lengths::empty>::create(window_dynamics{}, window_measurement{}, ctrlpp::nmhe_config<double, NX, NU, NY, window>{}));
        drive_past_warmup(estimator);

        REQUIRE(estimator.diagnostics().used_ekf_fallback);
        REQUIRE(estimator.diagnostics().status == ctrlpp::solve_status::invalid_backend_result);
    }
}

TEST_CASE("Nonlinear moving-horizon estimator accepts a conforming backend result unchanged", "[nmhe][result-shape][hardening]")
{
    auto estimator = ctrlpp::test::constructed(nonlinear_estimator<ctrlpp_test::report_lengths::conforming>::create(window_dynamics{}, window_measurement{}, ctrlpp::nmhe_config<double, NX, NU, NY, window>{}));
    drive_past_warmup(estimator);

    REQUIRE(!estimator.diagnostics().used_ekf_fallback);
    REQUIRE(estimator.diagnostics().status == ctrlpp::solve_status::optimal);
}
