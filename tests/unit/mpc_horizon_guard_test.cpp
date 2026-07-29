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

#include "ctrlpp/estimation/ekf.h"

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
#include <variant>

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

using constrained_nonlinear_estimator =
    ctrlpp::nmhe<double, NX, NU, NY, window,
                 ctrlpp_test::stub_nlp_solver<double>,
                 window_dynamics,
                 window_measurement,
                 1>;

template <ctrlpp_test::report_values Values, ctrlpp::solve_status Status>
using nonfinite_linear_estimator =
    ctrlpp::mhe<double, NX, NU, NY, window,
                ctrlpp_test::stub_qp_solver<double,
                                            ctrlpp_test::report_lengths::conforming,
                                            Values,
                                            Status>,
                window_dynamics,
                window_measurement>;

template <ctrlpp_test::report_values Values, ctrlpp::solve_status Status>
using nonfinite_nonlinear_estimator =
    ctrlpp::nmhe<double, NX, NU, NY, window,
                 ctrlpp_test::stub_nlp_solver<double,
                                              ctrlpp_test::report_lengths::conforming,
                                              Values,
                                              Status>,
                 window_dynamics,
                 window_measurement>;

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

template <ctrlpp_test::report_values Values, ctrlpp::solve_status Status>
void require_linear_estimator_rejects_nonfinite_result()
{
    auto estimator = ctrlpp::test::constructed(
        nonfinite_linear_estimator<Values, Status>::create(
            window_dynamics{},
            window_measurement{},
            ctrlpp::mhe_config<double, NX, NU, NY, window>{}));
    drive_past_warmup(estimator);
    REQUIRE(estimator.diagnostics().used_ekf_fallback);
    REQUIRE(estimator.diagnostics().status
            == ctrlpp::solve_status::invalid_backend_result);
    REQUIRE(estimator.state().allFinite());
}

template <ctrlpp_test::report_values Values, ctrlpp::solve_status Status>
void require_nonlinear_estimator_rejects_nonfinite_result()
{
    auto estimator = ctrlpp::test::constructed(
        nonfinite_nonlinear_estimator<Values, Status>::create(
            window_dynamics{},
            window_measurement{},
            ctrlpp::nmhe_config<double, NX, NU, NY, window>{}));
    drive_past_warmup(estimator);
    REQUIRE(estimator.diagnostics().used_ekf_fallback);
    REQUIRE(estimator.diagnostics().status
            == ctrlpp::solve_status::invalid_backend_result);
    REQUIRE(estimator.state().allFinite());
}

template <typename Estimator>
void require_refused_measurement_refills_window(Estimator& estimator)
{
    auto reference = ctrlpp::test::constructed(
        ctrlpp::ekf<double, NX, NU, NY,
                    window_dynamics,
                    window_measurement>::create(
            window_dynamics{},
            window_measurement{},
            ctrlpp::ekf_config<double, NX, NU, NY>{}));

    Eigen::Vector2d true_state = Eigen::Vector2d::Zero();
    auto advance = [&](std::size_t step) {
        Eigen::Matrix<double, 1, 1> input;
        auto const magnitude =
            static_cast<double>((step + 1) * (step + 1));
        input << ((step % 2 == 0) ? magnitude : -magnitude);
        true_state = window_dynamics{}(true_state, input);
        Eigen::Matrix<double, 1, 1> measurement;
        measurement << true_state(0);

        estimator.predict(input);
        reference.predict(input);
        REQUIRE(estimator.update(measurement).has_value());
        REQUIRE(reference.update(measurement).has_value());
    };

    for(std::size_t step = 0; step <= window; ++step)
        advance(step);
    REQUIRE(estimator.is_initialized());

    Eigen::Matrix<double, 1, 1> refused_input;
    refused_input << 100.0;
    true_state = window_dynamics{}(true_state, refused_input);
    estimator.predict(refused_input);
    reference.predict(refused_input);

    auto bad_measurement =
        Eigen::Matrix<double, 1, 1>::Constant(
            std::numeric_limits<double>::quiet_NaN());
    REQUIRE_FALSE(estimator.update(bad_measurement).has_value());
    REQUIRE_FALSE(reference.update(bad_measurement).has_value());
    REQUIRE_FALSE(estimator.is_initialized());
    REQUIRE(estimator.state() == reference.state());
    REQUIRE(estimator.covariance() == reference.covariance());

    for(std::size_t step = 0; step < window; ++step)
    {
        advance(window + 2 + step);
        REQUIRE(estimator.diagnostics().used_ekf_fallback);
        REQUIRE_FALSE(estimator.is_initialized());
        REQUIRE(estimator.state() == reference.state());
        REQUIRE(estimator.covariance() == reference.covariance());
    }

    advance(2 * window + 2);
    REQUIRE(estimator.is_initialized());
    REQUIRE_FALSE(estimator.diagnostics().used_ekf_fallback);
    REQUIRE(estimator.state().allFinite());
}

template <typename Result>
void require_moving_horizon_configuration_error(
    const Result& result,
    ctrlpp::moving_horizon_configuration_error expected)
{
    REQUIRE_FALSE(result.has_value());
    auto const* error =
        std::get_if<ctrlpp::moving_horizon_configuration_error>(&result.error());
    REQUIRE(error != nullptr);
    REQUIRE(*error == expected);
}

}

TEST_CASE("Refused moving-horizon measurements invalidate and refill the window",
          "[mhe][nmhe][window][hardening]")
{
    SECTION("linear estimator")
    {
        auto estimator = ctrlpp::test::constructed(
            linear_estimator<ctrlpp_test::report_lengths::conforming>::create(
                window_dynamics{},
                window_measurement{},
                ctrlpp::mhe_config<double, NX, NU, NY, window>{}));
        require_refused_measurement_refills_window(estimator);
    }

    SECTION("nonlinear estimator")
    {
        auto estimator = ctrlpp::test::constructed(
            nonlinear_estimator<ctrlpp_test::report_lengths::conforming>::create(
                window_dynamics{},
                window_measurement{},
                ctrlpp::nmhe_config<double, NX, NU, NY, window>{}));
        require_refused_measurement_refills_window(estimator);
    }
}

TEST_CASE("Moving-horizon factories reject singular formulation weights",
          "[mhe][nmhe][configuration][hardening]")
{
    auto require_both = [](auto configure,
                           ctrlpp::moving_horizon_configuration_error error) {
        ctrlpp::mhe_config<double, NX, NU, NY, window> linear;
        configure(linear);
        require_moving_horizon_configuration_error(
            linear_estimator<ctrlpp_test::report_lengths::conforming>::create(
                window_dynamics{}, window_measurement{}, linear),
            error);

        ctrlpp::nmhe_config<double, NX, NU, NY, window> nonlinear;
        configure(nonlinear);
        require_moving_horizon_configuration_error(
            nonlinear_estimator<ctrlpp_test::report_lengths::conforming>::create(
                window_dynamics{}, window_measurement{}, nonlinear),
            error);
    };

    require_both(
        [](auto& config) { config.Q.setZero(); },
        ctrlpp::moving_horizon_configuration_error::
            non_invertible_process_noise);
    require_both(
        [](auto& config) { config.R.setZero(); },
        ctrlpp::moving_horizon_configuration_error::
            non_invertible_measurement_noise);
    require_both(
        [](auto& config) { config.P0.setZero(); },
        ctrlpp::moving_horizon_configuration_error::
            non_invertible_initial_covariance);
}

TEST_CASE("Moving-horizon factories preserve embedded-filter configuration errors",
          "[mhe][nmhe][configuration][hardening]")
{
    auto require_filter_error = [](const auto& result) {
        REQUIRE_FALSE(result.has_value());
        auto const* error = std::get_if<ctrlpp::filter_error>(&result.error());
        REQUIRE(error != nullptr);
        REQUIRE(*error == ctrlpp::filter_error::non_finite_process_noise);
    };

    ctrlpp::mhe_config<double, NX, NU, NY, window> linear;
    linear.Q(0, 0) = std::numeric_limits<double>::quiet_NaN();
    require_filter_error(
        linear_estimator<ctrlpp_test::report_lengths::conforming>::create(
            window_dynamics{}, window_measurement{}, linear));

    ctrlpp::nmhe_config<double, NX, NU, NY, window> nonlinear;
    nonlinear.Q(0, 0) = std::numeric_limits<double>::quiet_NaN();
    require_filter_error(
        nonlinear_estimator<ctrlpp_test::report_lengths::conforming>::create(
            window_dynamics{}, window_measurement{}, nonlinear));
}

TEST_CASE("Moving-horizon factories validate every common formulation option",
          "[mhe][nmhe][configuration][hardening]")
{
    auto require_both = [](auto configure,
                           ctrlpp::moving_horizon_configuration_error error) {
        ctrlpp::mhe_config<double, NX, NU, NY, window> linear;
        configure(linear);
        require_moving_horizon_configuration_error(
            linear_estimator<ctrlpp_test::report_lengths::conforming>::create(
                window_dynamics{}, window_measurement{}, linear),
            error);

        ctrlpp::nmhe_config<double, NX, NU, NY, window> nonlinear;
        configure(nonlinear);
        require_moving_horizon_configuration_error(
            nonlinear_estimator<ctrlpp_test::report_lengths::conforming>::create(
                window_dynamics{}, window_measurement{}, nonlinear),
            error);
    };

    auto const infinity = std::numeric_limits<double>::infinity();
    auto const nan = std::numeric_limits<double>::quiet_NaN();

    require_both(
        [](auto& config) { config.arrival_cost_weight = 0.0; },
        ctrlpp::moving_horizon_configuration_error::
            non_positive_arrival_cost_weight);
    require_both(
        [infinity](auto& config) {
            config.arrival_cost_weight = infinity;
        },
        ctrlpp::moving_horizon_configuration_error::
            non_positive_arrival_cost_weight);
    require_both(
        [infinity](auto& config) {
            config.x_min = config.x0;
            config.x_min->setConstant(infinity);
        },
        ctrlpp::moving_horizon_configuration_error::invalid_state_bounds);
    require_both(
        [](auto& config) {
            config.x_min = config.x0;
            config.x_max = config.x0;
            (*config.x_min)(0) = 1.0;
            (*config.x_max)(0) = -1.0;
        },
        ctrlpp::moving_horizon_configuration_error::invalid_state_bounds);
    require_both(
        [](auto& config) {
            config.residual_bound = config.R.col(0);
            config.residual_bound->setConstant(-1.0);
        },
        ctrlpp::moving_horizon_configuration_error::invalid_residual_bound);
    require_both(
        [infinity](auto& config) {
            config.residual_bound = config.R.col(0);
            config.residual_bound->setConstant(infinity);
        },
        ctrlpp::moving_horizon_configuration_error::invalid_residual_bound);
    require_both(
        [](auto& config) { config.soft_penalty = 0.0; },
        ctrlpp::moving_horizon_configuration_error::
            non_positive_soft_penalty);
    require_both(
        [infinity](auto& config) { config.soft_penalty = infinity; },
        ctrlpp::moving_horizon_configuration_error::
            non_positive_soft_penalty);
    require_both(
        [](auto& config) { config.numerical_eps = 0.0; },
        ctrlpp::moving_horizon_configuration_error::
            non_positive_numerical_eps);
    require_both(
        [nan](auto& config) { config.numerical_eps = nan; },
        ctrlpp::moving_horizon_configuration_error::
            non_positive_numerical_eps);
}

TEST_CASE("Nonlinear moving-horizon factory validates path-constraint operands",
          "[nmhe][configuration][hardening]")
{
    ctrlpp::nmhe_config<double, NX, NU, NY, window, 1> config;
    config.path_constraint = [](const Eigen::Vector2d&) {
        return Eigen::Matrix<double, 1, 1>::Zero();
    };

    config.path_penalty << 0.0;
    require_moving_horizon_configuration_error(
        constrained_nonlinear_estimator::create(
            window_dynamics{}, window_measurement{}, config),
        ctrlpp::moving_horizon_configuration_error::invalid_path_penalty);

    config.path_penalty << std::numeric_limits<double>::infinity();
    require_moving_horizon_configuration_error(
        constrained_nonlinear_estimator::create(
            window_dynamics{}, window_measurement{}, config),
        ctrlpp::moving_horizon_configuration_error::invalid_path_penalty);

    config.path_penalty << 1.0;
    config.path_constraint = [](const Eigen::Vector2d&) {
        return Eigen::Matrix<double, 1, 1>::Constant(
            std::numeric_limits<double>::quiet_NaN());
    };
    require_moving_horizon_configuration_error(
        constrained_nonlinear_estimator::create(
            window_dynamics{}, window_measurement{}, config),
        ctrlpp::moving_horizon_configuration_error::
            non_finite_path_constraint);
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

TEST_CASE("Linear MPC rejects non-finite results under every accepted status",
          "[mpc][result-finiteness][hardening]")
{
    for(auto const status : {ctrlpp::solve_status::optimal,
                             ctrlpp::solve_status::solved_inaccurate,
                             ctrlpp::solve_status::max_iterations,
                             ctrlpp::solve_status::time_limit})
    {
        for(auto const values : {ctrlpp_test::report_values::nan,
                                 ctrlpp_test::report_values::positive_infinity,
                                 ctrlpp_test::report_values::negative_infinity})
        {
            auto solver = ctrlpp_test::stub_qp_solver<double>{
                ctrlpp_test::report_lengths::conforming, values, status};
            auto controller =
                linear_controller::create(make_system(),
                                          linear_config(shape_horizon),
                                          std::move(solver));
            REQUIRE(controller.has_value());

            auto solved = controller->solve(Eigen::Vector2d{1.0, 0.0});
            REQUIRE_FALSE(solved.has_value());
            REQUIRE(solved.error()
                    == ctrlpp::solver_error::invalid_backend_result);
            REQUIRE(controller->diagnostics().status
                    == ctrlpp::solve_status::invalid_backend_result);
            REQUIRE_FALSE(controller->trajectory().has_value());
        }
    }
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

TEST_CASE("Linear moving-horizon estimator rejects non-finite accepted results",
          "[mhe][result-finiteness][hardening]")
{
    require_linear_estimator_rejects_nonfinite_result<
        ctrlpp_test::report_values::nan,
        ctrlpp::solve_status::optimal>();
    require_linear_estimator_rejects_nonfinite_result<
        ctrlpp_test::report_values::positive_infinity,
        ctrlpp::solve_status::optimal>();
    require_linear_estimator_rejects_nonfinite_result<
        ctrlpp_test::report_values::negative_infinity,
        ctrlpp::solve_status::optimal>();
    require_linear_estimator_rejects_nonfinite_result<
        ctrlpp_test::report_values::nan,
        ctrlpp::solve_status::solved_inaccurate>();
    require_linear_estimator_rejects_nonfinite_result<
        ctrlpp_test::report_values::positive_infinity,
        ctrlpp::solve_status::solved_inaccurate>();
    require_linear_estimator_rejects_nonfinite_result<
        ctrlpp_test::report_values::negative_infinity,
        ctrlpp::solve_status::solved_inaccurate>();
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

TEST_CASE("Nonlinear moving-horizon estimator rejects non-finite accepted results",
          "[nmhe][result-finiteness][hardening]")
{
    require_nonlinear_estimator_rejects_nonfinite_result<
        ctrlpp_test::report_values::nan,
        ctrlpp::solve_status::optimal>();
    require_nonlinear_estimator_rejects_nonfinite_result<
        ctrlpp_test::report_values::positive_infinity,
        ctrlpp::solve_status::optimal>();
    require_nonlinear_estimator_rejects_nonfinite_result<
        ctrlpp_test::report_values::negative_infinity,
        ctrlpp::solve_status::optimal>();
    require_nonlinear_estimator_rejects_nonfinite_result<
        ctrlpp_test::report_values::nan,
        ctrlpp::solve_status::solved_inaccurate>();
    require_nonlinear_estimator_rejects_nonfinite_result<
        ctrlpp_test::report_values::positive_infinity,
        ctrlpp::solve_status::solved_inaccurate>();
    require_nonlinear_estimator_rejects_nonfinite_result<
        ctrlpp_test::report_values::negative_infinity,
        ctrlpp::solve_status::solved_inaccurate>();
}
