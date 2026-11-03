// Guards on the same data path: the horizon that sizes the posed problem, the
// configuration that shapes it, and the shape of the result the backend hands
// back.
//
// Construction-time horizon domain. The prediction horizon is the one runtime
// quantity that scales every derived dimension of the posed problem, and it is
// caller-supplied. These cases pin the boundary of its accepted domain from both
// sides and assert the specific typed error on each rejection, never merely that
// construction "did not crash".
//
// Configuration domain. Bounds, the rate bound, the cost weights and the
// constraint penalties are each pinned from both sides of their boundary, on
// the linear and the runtime-horizon nonlinear factory, through both overloads
// of each. The linear factory's plant and the weights it maps through the
// plant's output matrix are pinned the same way. The compile-time-horizon factory shares the admission; its cases pin
// the order of admission, formulation refusal and solver setup.
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
#include "scripted_nlp_solver.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>

#include <span>
#include <array>
#include <cmath>
#include <tuple>
#include <limits>
#include <memory>
#include <random>
#include <vector>
#include <cstddef>
#include <cstring>
#include <utility>
#include <variant>
#include <optional>
#include <algorithm>
#include <functional>
#include <type_traits>

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

// The Jacobians are exact, so the fallback filter's linearization carries no
// finite-difference step whose rounding it would amplify.
struct window_dynamics
{
    auto operator()(const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>& u) const -> Eigen::Vector2d { return Eigen::Vector2d{x(0) + dt * x(1), x(1) + dt * u(0)}; }

    auto jacobian_x(const Eigen::Vector2d&, const Eigen::Matrix<double, 1, 1>&) const -> Eigen::Matrix2d
    {
        Eigen::Matrix2d F;
        F << 1.0, dt, 0.0, 1.0;
        return F;
    }

    auto jacobian_u(const Eigen::Vector2d&, const Eigen::Matrix<double, 1, 1>&) const -> Eigen::Vector2d { return Eigen::Vector2d{0.0, dt}; }
};

struct window_measurement
{
    auto operator()(const Eigen::Vector2d& x) const -> Eigen::Matrix<double, 1, 1> { return x.head<1>(); }

    auto jacobian(const Eigen::Vector2d&) const -> Eigen::Matrix<double, 1, 2> { return Eigen::Matrix<double, 1, 2>{1.0, 0.0}; }
};

// How far two runs of the fallback filter on the same inputs may part, to first
// order in eps. They part only by rounding, which a compiler may perform
// differently for each, fusing multiply-adds in one inlining context and not in
// another. Each step's partings are propagated through the reference run's own
// recursion and the step's roundings are added at eps, which covers both runs:
// the state prediction rounds each component twice, the covariance prediction
// six times (two two-term inner products, the noise sum and the
// symmetrization), the gain divides once by an S rounded once, the correction
// rounds three times, and the Joseph update seven (forming I - K H, two two-term
// inner products, the noise sum and the symmetrization). The Joseph form is
// first-order insensitive to the gain, so a gain parting moves the state but
// not the covariance.
struct fallback_parting
{
    Eigen::Vector2d state{Eigen::Vector2d::Zero()};
    Eigen::Matrix2d covariance{Eigen::Matrix2d::Zero()};

    void predict(const Eigen::Vector2d& x, const Eigen::Matrix2d& P, const Eigen::Matrix2d& Q, const Eigen::Matrix<double, 1, 1>& u)
    {
        const double eps = std::numeric_limits<double>::epsilon();
        const Eigen::Matrix2d F = window_dynamics{}.jacobian_x(x, u).cwiseAbs();
        const Eigen::Vector2d G = window_dynamics{}.jacobian_u(x, u).cwiseAbs();
        state = (F * state + 2.0 * eps * (F * x.cwiseAbs() + G * std::abs(u(0)))).eval();
        covariance = (F * covariance * F.transpose() + 6.0 * eps * (F * P.cwiseAbs() * F.transpose() + Q.cwiseAbs())).eval();
    }

    void update(const Eigen::Vector2d& x_prior, const Eigen::Matrix2d& P_prior, double innovation, double r)
    {
        const double eps = std::numeric_limits<double>::epsilon();
        const double S = P_prior(0, 0) + r;
        const Eigen::Vector2d K = P_prior.col(0) / S;
        const Eigen::Vector2d K_parting = (covariance.col(0) + K.cwiseAbs() * covariance(0, 0)) / S + 2.0 * eps * K.cwiseAbs();
        const double innovation_parting = state(0) + eps * std::abs(innovation);
        state = (state + K_parting * std::abs(innovation) + K.cwiseAbs() * innovation_parting
                 + 3.0 * eps * (x_prior.cwiseAbs() + K.cwiseAbs() * std::abs(innovation))).eval();
        Eigen::Matrix2d IKH = Eigen::Matrix2d::Identity();
        IKH.col(0) -= K;
        const Eigen::Matrix2d A = IKH.cwiseAbs();
        covariance = (A * covariance * A.transpose()
                      + 7.0 * eps * (A * P_prior.cwiseAbs() * A.transpose() + r * K.cwiseAbs() * K.cwiseAbs().transpose())).eval();
    }
};

template <typename Estimator, typename Reference>
void require_within_parting(const Estimator& estimator, const Reference& reference, const fallback_parting& parting)
{
    CAPTURE(estimator.state(), reference.state(), parting.state, estimator.covariance(), reference.covariance(), parting.covariance);
    REQUIRE(((estimator.state() - reference.state()).cwiseAbs().array() <= parting.state.array()).all());
    REQUIRE(((estimator.covariance() - reference.covariance()).cwiseAbs().array() <= parting.covariance.array()).all());
}

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

    const ctrlpp::ekf_config<double, NX, NU, NY> noise{};
    const double r = noise.R(0, 0);
    fallback_parting parting;
    auto predict_both = [&](const Eigen::Matrix<double, 1, 1>& input) {
        const Eigen::Vector2d x = reference.state();
        const Eigen::Matrix2d P = reference.covariance();
        estimator.predict(input);
        reference.predict(input);
        parting.predict(x, P, noise.Q, input);
    };

    Eigen::Vector2d true_state = Eigen::Vector2d::Zero();
    auto advance = [&](std::size_t step) {
        Eigen::Matrix<double, 1, 1> input;
        auto const magnitude =
            static_cast<double>((step + 1) * (step + 1));
        input << ((step % 2 == 0) ? magnitude : -magnitude);
        true_state = window_dynamics{}(true_state, input);
        Eigen::Matrix<double, 1, 1> measurement;
        measurement << true_state(0);

        predict_both(input);
        const Eigen::Vector2d x_prior = reference.state();
        const Eigen::Matrix2d P_prior = reference.covariance();
        REQUIRE(estimator.update(measurement).has_value());
        REQUIRE(reference.update(measurement).has_value());
        parting.update(x_prior, P_prior, measurement(0) - x_prior(0), r);
    };

    for(std::size_t step = 0; step <= window; ++step)
        advance(step);
    REQUIRE(estimator.is_initialized());

    Eigen::Matrix<double, 1, 1> refused_input;
    refused_input << 100.0;
    true_state = window_dynamics{}(true_state, refused_input);
    predict_both(refused_input);

    auto bad_measurement =
        Eigen::Matrix<double, 1, 1>::Constant(
            std::numeric_limits<double>::quiet_NaN());
    REQUIRE_FALSE(estimator.update(bad_measurement).has_value());
    REQUIRE_FALSE(reference.update(bad_measurement).has_value());
    REQUIRE_FALSE(estimator.is_initialized());
    require_within_parting(estimator, reference, parting);

    for(std::size_t step = 0; step < window; ++step)
    {
        advance(window + 2 + step);
        REQUIRE(estimator.diagnostics().used_ekf_fallback);
        REQUIRE_FALSE(estimator.is_initialized());
        require_within_parting(estimator, reference, parting);
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

template <typename Alternative, typename Result>
void require_construction_error(const Result& result, Alternative expected)
{
    REQUIRE_FALSE(result.has_value());
    auto const* error = std::get_if<Alternative>(&result.error());
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

TEST_CASE("Moving-horizon factories report an indefinite process noise as the embedded filter's error",
          "[mhe][nmhe][configuration][hardening]")
{
    auto require_indefinite = [](const auto& result) {
        REQUIRE_FALSE(result.has_value());
        auto const* error = std::get_if<ctrlpp::filter_error>(&result.error());
        REQUIRE(error != nullptr);
        REQUIRE(*error == ctrlpp::filter_error::indefinite_process_noise);
    };

    ctrlpp::mhe_config<double, NX, NU, NY, window> linear;
    linear.Q = Eigen::Vector2d{1.0, -1.0}.asDiagonal();
    require_indefinite(
        linear_estimator<ctrlpp_test::report_lengths::conforming>::create(
            window_dynamics{}, window_measurement{}, linear));

    ctrlpp::nmhe_config<double, NX, NU, NY, window> nonlinear;
    nonlinear.Q = Eigen::Vector2d{1.0, -1.0}.asDiagonal();
    require_indefinite(
        nonlinear_estimator<ctrlpp_test::report_lengths::conforming>::create(
            window_dynamics{}, window_measurement{}, nonlinear));
}

TEST_CASE("Moving-horizon factories report a structurally asymmetric process noise as the embedded filter's error",
          "[mhe][nmhe][configuration][hardening]")
{
    auto require_asymmetric = [](const auto& result) {
        REQUIRE_FALSE(result.has_value());
        auto const* error = std::get_if<ctrlpp::filter_error>(&result.error());
        REQUIRE(error != nullptr);
        REQUIRE(*error == ctrlpp::filter_error::asymmetric_process_noise);
    };

    ctrlpp::mhe_config<double, NX, NU, NY, window> linear;
    linear.Q << 1.0, 0.5, 0.4, 1.0;
    require_asymmetric(
        linear_estimator<ctrlpp_test::report_lengths::conforming>::create(
            window_dynamics{}, window_measurement{}, linear));

    ctrlpp::nmhe_config<double, NX, NU, NY, window> nonlinear;
    nonlinear.Q << 1.0, 0.5, 0.4, 1.0;
    require_asymmetric(
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

    require_construction_error(result, ctrlpp::controller_construction_error::non_positive_horizon);
}

TEST_CASE("Linear MPC rejects a zero horizon at construction", "[mpc][horizon][hardening]")
{
    auto result = linear_controller::create(make_system(), linear_config(0), ctrlpp_test::stub_qp_solver<double>{});

    require_construction_error(result, ctrlpp::controller_construction_error::non_positive_horizon);
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
    require_construction_error(just_past, ctrlpp::controller_construction_error::horizon_overflow);

    auto largest = linear_controller::create(make_system(), linear_config(std::numeric_limits<int>::max()), ctrlpp_test::stub_qp_solver<double>{});
    require_construction_error(largest, ctrlpp::controller_construction_error::horizon_overflow);
}

TEST_CASE("Linear MPC reports a solver setup failure at construction with the backend's cause", "[mpc][setup][hardening]")
{
    ctrlpp_test::stub_qp_solver<double> solver;
    solver.fail_setup = true;
    auto result = linear_controller::create(make_system(), linear_config(shape_horizon), std::move(solver));

    require_construction_error(result, ctrlpp_test::stub_setup_error::rejected);
}

TEST_CASE("Nonlinear MPC rejects a negative horizon at construction", "[nmpc][horizon][hardening]")
{
    auto result = nonlinear_controller::create(double_integrator, nonlinear_config(-1), ctrlpp_test::stub_nlp_solver<double>{});
    require_construction_error(result, ctrlpp::controller_construction_error::non_positive_horizon);
}

TEST_CASE("Nonlinear MPC rejects a zero horizon at construction", "[nmpc][horizon][hardening]")
{
    auto result = nonlinear_controller::create(double_integrator, nonlinear_config(0), ctrlpp_test::stub_nlp_solver<double>{});
    require_construction_error(result, ctrlpp::controller_construction_error::non_positive_horizon);
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
    require_construction_error(just_past, ctrlpp::controller_construction_error::horizon_overflow);

    auto largest = nonlinear_controller::create(double_integrator, nonlinear_config(std::numeric_limits<int>::max()), ctrlpp_test::stub_nlp_solver<double>{});
    require_construction_error(largest, ctrlpp::controller_construction_error::horizon_overflow);
}

TEST_CASE("Runtime-horizon nonlinear MPC reports a solver setup failure at construction with the backend's cause", "[nmpc][setup][hardening]")
{
    ctrlpp_test::stub_nlp_solver<double> solver;
    solver.fail_setup = true;
    auto result = nonlinear_controller::create(double_integrator, nonlinear_config(shape_horizon), std::move(solver));

    require_construction_error(result, ctrlpp_test::stub_setup_error::rejected);
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

TEST_CASE("Compile-time-horizon controller refuses a disagreeing runtime horizon at construction", "[nmpc][horizon][hardening]")
{
    constexpr std::size_t NH = 4;
    using static_controller = ctrlpp::nmpc_static<double, NX, NU, NH, ctrlpp_test::stub_nlp_solver<double>, decltype(double_integrator)>;

    for(const int runtime_horizon : {static_cast<int>(NH) - 1, static_cast<int>(NH) + 1, 0})
    {
        CAPTURE(runtime_horizon);
        require_construction_error(static_controller::create(double_integrator, nonlinear_config(runtime_horizon)), ctrlpp::nlp_formulation_error::horizon_mismatch);
        require_construction_error(static_controller::create(double_integrator, nonlinear_config(runtime_horizon), ctrlpp_test::stub_nlp_solver<double>{}),
                                             ctrlpp::nlp_formulation_error::horizon_mismatch);
    }
}

namespace
{

using constrained_default_alias = ctrlpp::nmpc<double, NX, NU, 5, ctrlpp_test::stub_nlp_solver<double>, decltype(double_integrator), 1>;

auto path_constrained_config() -> ctrlpp::nmpc_config<double, NX, NU, 1>
{
    ctrlpp::nmpc_config<double, NX, NU, 1> cfg;
    cfg.horizon = shape_horizon;
    cfg.path_constraint = [](const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>&) { return Eigen::Matrix<double, 1, 1>{x(0) - 10.0}; };
    return cfg;
}

}

TEST_CASE("The default nmpc refuses the default configuration with a path constraint at construction", "[nmpc][configuration][hardening]")
{
    const auto cfg = path_constrained_config();
    REQUIRE(cfg.soft_constraints);

    require_construction_error(constrained_default_alias::create(double_integrator, cfg), ctrlpp::nlp_formulation_error::slack_not_supported);
    require_construction_error(constrained_default_alias::create(double_integrator, cfg, ctrlpp_test::stub_nlp_solver<double>{}),
                                         ctrlpp::nlp_formulation_error::slack_not_supported);
}

TEST_CASE("The default nmpc solves the configurations it can pose", "[nmpc][configuration][hardening]")
{
    SECTION("a path constraint with soft constraints turned off")
    {
        auto cfg = path_constrained_config();
        cfg.soft_constraints = false;
        auto controller = constrained_default_alias::create(double_integrator, cfg);
        REQUIRE(controller.has_value());
        REQUIRE(controller->solve(Eigen::Vector2d{1.0, 0.0}).has_value());
    }

    SECTION("the default configuration without a constraint")
    {
        ctrlpp::nmpc_config<double, NX, NU, 1> cfg;
        cfg.horizon = shape_horizon;
        REQUIRE(cfg.soft_constraints);
        auto controller = constrained_default_alias::create(double_integrator, cfg, ctrlpp_test::stub_nlp_solver<double>{});
        REQUIRE(controller.has_value());
        REQUIRE(controller->solve(Eigen::Vector2d{1.0, 0.0}).has_value());
    }
}

TEST_CASE("Compile-time-horizon nonlinear MPC reports a solver setup failure at construction with the backend's cause", "[nmpc][setup][hardening]")
{
    ctrlpp_test::stub_nlp_solver<double> solver;
    solver.fail_setup = true;
    auto cfg = path_constrained_config();
    cfg.soft_constraints = false;
    require_construction_error(constrained_default_alias::create(double_integrator, cfg, std::move(solver)), ctrlpp_test::stub_setup_error::rejected);
}

TEST_CASE("Compile-time-horizon nonlinear MPC checks admission, then the formulation, then setup", "[nmpc][configuration][hardening]")
{
    auto cfg = path_constrained_config();
    cfg.horizon = shape_horizon + 1;

    SECTION("an inverted input bound is reported before the formulation refusals")
    {
        cfg.u_min = Eigen::Matrix<double, 1, 1>{1.0};
        cfg.u_max = Eigen::Matrix<double, 1, 1>{-1.0};
        ctrlpp_test::stub_nlp_solver<double> refusing;
        refusing.fail_setup = true;
        require_construction_error(constrained_default_alias::create(double_integrator, cfg), ctrlpp::controller_construction_error::invalid_input_bounds);
        require_construction_error(constrained_default_alias::create(double_integrator, cfg, std::move(refusing)),
                                             ctrlpp::controller_construction_error::invalid_input_bounds);
    }

    SECTION("a horizon mismatch is reported before the slack refusal")
    {
        require_construction_error(constrained_default_alias::create(double_integrator, cfg), ctrlpp::nlp_formulation_error::horizon_mismatch);
    }

    SECTION("a formulation refusal is reported before the solver is set up")
    {
        cfg.horizon = shape_horizon;
        ctrlpp_test::stub_nlp_solver<double> refusing;
        refusing.fail_setup = true;
        require_construction_error(constrained_default_alias::create(double_integrator, cfg, std::move(refusing)),
                                             ctrlpp::nlp_formulation_error::slack_not_supported);
    }
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
        auto controller = ctrlpp::test::constructed(compile_time_controller<ctrlpp_test::report_lengths::short_primal>::create(double_integrator, nonlinear_config(shape_horizon)));

        auto solved = controller.solve(Eigen::Vector2d{1.0, 0.0});
        REQUIRE(!solved.has_value());
        REQUIRE(solved.error() == ctrlpp::solver_error::invalid_backend_result);
    }

    SECTION("a primal that is empty despite the reported optimal status")
    {
        auto controller = ctrlpp::test::constructed(compile_time_controller<ctrlpp_test::report_lengths::empty>::create(double_integrator, nonlinear_config(shape_horizon)));

        auto solved = controller.solve(Eigen::Vector2d{1.0, 0.0});
        REQUIRE(!solved.has_value());
        REQUIRE(solved.error() == ctrlpp::solver_error::invalid_backend_result);
    }
}

TEST_CASE("Compile-time-horizon nonlinear MPC accepts a conforming backend result unchanged", "[nmpc][result-shape][hardening]")
{
    auto controller = ctrlpp::test::constructed(compile_time_controller<ctrlpp_test::report_lengths::conforming>::create(double_integrator, nonlinear_config(shape_horizon)));

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

namespace
{

template <bool SolveInto>
using scripted_static_controller = ctrlpp::nmpc_static<double, NX, NU, 5, ctrlpp_test::scripted_nlp_solver<double, SolveInto>, decltype(double_integrator)>;

using scripted_dynamic_controller = ctrlpp::nmpc_dynamic<double, NX, NU, ctrlpp_test::scripted_nlp_solver<double>, decltype(double_integrator)>;

using scripted_record = ctrlpp_test::scripted_record<double>;

constexpr std::array accepted_statuses{ctrlpp::solve_status::optimal, ctrlpp::solve_status::solved_inaccurate, ctrlpp::solve_status::max_iterations, ctrlpp::solve_status::time_limit};

struct corruption
{
    ctrlpp_test::corrupted_field field;
    ctrlpp_test::report_values value;
};

// Every value-carrying field alone with each non-finite value, then the iteration
// count alone, which the stub sets to -1 whatever the value.
auto single_field_corruptions() -> std::vector<corruption>
{
    std::vector<corruption> out;
    for(auto const field : {ctrlpp_test::corrupted_field::primal, ctrlpp_test::corrupted_field::objective, ctrlpp_test::corrupted_field::solve_time, ctrlpp_test::corrupted_field::primal_residual})
        for(auto const value : {ctrlpp_test::report_values::nan, ctrlpp_test::report_values::positive_infinity, ctrlpp_test::report_values::negative_infinity})
            out.push_back({field, value});
    out.push_back({ctrlpp_test::corrupted_field::iterations, ctrlpp_test::report_values::nan});
    return out;
}

auto scripted(std::vector<ctrlpp_test::scripted_answer> script) -> std::shared_ptr<scripted_record>
{
    auto record = std::make_shared<scripted_record>();
    record->script = std::move(script);
    return record;
}

auto rate_limited_config() -> ctrlpp::nmpc_config<double, NX, NU>
{
    auto cfg = nonlinear_config(shape_horizon);
    cfg.du_max = Eigen::Matrix<double, 1, 1>{0.5};
    return cfg;
}

auto make_scripted_dynamic(std::shared_ptr<scripted_record> record, const ctrlpp::nmpc_config<double, NX, NU>& cfg) -> scripted_dynamic_controller
{
    return ctrlpp::test::constructed(scripted_dynamic_controller::create(double_integrator, cfg, ctrlpp_test::scripted_nlp_solver<double>{std::move(record)}));
}

template <bool SolveInto>
auto make_scripted_static(std::shared_ptr<scripted_record> record, const ctrlpp::nmpc_config<double, NX, NU>& cfg) -> scripted_static_controller<SolveInto>
{
    return ctrlpp::test::constructed(scripted_static_controller<SolveInto>::create(double_integrator, cfg, ctrlpp_test::scripted_nlp_solver<double, SolveInto>{std::move(record)}));
}

auto bitwise_equal(const Eigen::VectorXd& a, const Eigen::VectorXd& b) -> bool
{
    return a.size() == b.size() && std::memcmp(a.data(), b.data(), sizeof(double) * static_cast<std::size_t>(a.size())) == 0;
}

template <typename Controller>
auto constraint_rows_at_origin(const Controller& controller) -> Eigen::VectorXd
{
    const auto& problem = controller.problem();
    const Eigen::VectorXd z = Eigen::VectorXd::Zero(problem.n_vars);
    Eigen::VectorXd rows = Eigen::VectorXd::Zero(problem.n_constraints);
    problem.constraints_cb(std::span<const double>{z.data(), static_cast<std::size_t>(z.size())}, std::span<double>{rows.data(), static_cast<std::size_t>(rows.size())});
    return rows;
}

template <typename Controller>
void check_refused(Controller& controller)
{
    auto solved = controller.solve(Eigen::Vector2d{1.0, 0.0});
    CHECK_FALSE(solved.has_value());
    if(!solved.has_value())
        CHECK(solved.error() == ctrlpp::solver_error::invalid_backend_result);
    CHECK(controller.diagnostics().status == ctrlpp::solve_status::invalid_backend_result);
    if constexpr(requires { controller.trajectory(); })
        CHECK_FALSE(controller.trajectory().has_value());
}

template <typename Make>
void check_every_field_refused(Make make)
{
    for(auto const status : accepted_statuses)
    {
        for(auto const bad : single_field_corruptions())
        {
            CAPTURE(status, bad.field, bad.value);
            auto controller = make(scripted({ctrlpp_test::scripted_answer{status, bad.field, bad.value}}));
            check_refused(controller);
        }
    }
}

// Script: finite, corrupted, finite. The refused middle call must leave the next
// call's warm start, the previous input inside the posed rate rows, and the
// stored primal exactly as the first accepted call left them. The rows are
// compared after the second and the third call, because each solve poses the
// previous input before it calls the backend; the first call's rows guard that
// the comparison is not vacuous.
template <typename Make>
void check_refusal_commits_nothing(Make make)
{
    const Eigen::Vector2d x0{1.0, 0.0};
    for(auto const status : accepted_statuses)
    {
        for(auto const bad : single_field_corruptions())
        {
            CAPTURE(status, bad.field, bad.value);
            const ctrlpp_test::scripted_answer finite{status};
            auto record = scripted({finite, ctrlpp_test::scripted_answer{status, bad.field, bad.value}, finite});
            auto controller = make(record);

            REQUIRE(controller.solve(x0).has_value());
            const Eigen::VectorXd rows_first = constraint_rows_at_origin(controller);
            const Eigen::VectorXd primal_first = controller.last_solution();
            CHECK(bitwise_equal(primal_first, record->primals[0]));

            CHECK_FALSE(controller.solve(x0).has_value());
            const Eigen::VectorXd rows_refused = constraint_rows_at_origin(controller);
            CHECK(bitwise_equal(controller.last_solution(), primal_first));

            CHECK(controller.solve(x0).has_value());
            REQUIRE(record->warm_starts.size() == 3);
            CHECK(bitwise_equal(record->warm_starts[2], record->warm_starts[1]));
            CHECK(rows_refused != rows_first);
            CHECK(bitwise_equal(constraint_rows_at_origin(controller), rows_refused));
        }
    }
}

}

TEST_CASE("Runtime-horizon nonlinear MPC refuses a NaN primal under every accepted status", "[nmpc][result-finiteness][hardening]")
{
    for(auto const status : accepted_statuses)
    {
        CAPTURE(status);
        auto controller = make_scripted_dynamic(scripted({ctrlpp_test::scripted_answer{status, ctrlpp_test::corrupted_field::primal}}), nonlinear_config(shape_horizon));

        auto solved = controller.solve(Eigen::Vector2d{1.0, 0.0});
        REQUIRE_FALSE(solved.has_value());
        REQUIRE(solved.error() == ctrlpp::solver_error::invalid_backend_result);
        REQUIRE(controller.diagnostics().status == ctrlpp::solve_status::invalid_backend_result);
        REQUIRE_FALSE(controller.trajectory().has_value());
    }
}

TEST_CASE("Runtime-horizon nonlinear MPC refuses each malformed result field under every accepted status", "[nmpc][result-finiteness][hardening]")
{
    check_every_field_refused([](std::shared_ptr<scripted_record> record) { return make_scripted_dynamic(std::move(record), nonlinear_config(shape_horizon)); });
}

TEST_CASE("Compile-time-horizon nonlinear MPC refuses each malformed result field under every accepted status", "[nmpc][result-finiteness][hardening]")
{
    SECTION("a backend that returns its primal by value")
    {
        check_every_field_refused([](std::shared_ptr<scripted_record> record) { return make_scripted_static<false>(std::move(record), nonlinear_config(shape_horizon)); });
    }

    SECTION("a backend that writes its primal into the controller's buffer")
    {
        check_every_field_refused([](std::shared_ptr<scripted_record> record) { return make_scripted_static<true>(std::move(record), nonlinear_config(shape_horizon)); });
    }
}

TEST_CASE("Runtime-horizon nonlinear MPC commits nothing from a refused result", "[nmpc][result-finiteness][hardening]")
{
    check_refusal_commits_nothing([](std::shared_ptr<scripted_record> record) { return make_scripted_dynamic(std::move(record), rate_limited_config()); });
}

TEST_CASE("Compile-time-horizon nonlinear MPC commits nothing from a refused result", "[nmpc][result-finiteness][hardening]")
{
    SECTION("a backend that returns its primal by value")
    {
        check_refusal_commits_nothing([](std::shared_ptr<scripted_record> record) { return make_scripted_static<false>(std::move(record), rate_limited_config()); });
    }

    SECTION("a backend that writes its primal into the controller's buffer")
    {
        check_refusal_commits_nothing([](std::shared_ptr<scripted_record> record) { return make_scripted_static<true>(std::move(record), rate_limited_config()); });
    }
}

TEST_CASE("Compile-time-horizon nonlinear MPC keeps its stored primal at the posed dimension on a by-value backend", "[nmpc][result-finiteness][hardening]")
{
    constexpr auto posed = static_cast<Eigen::Index>(scripted_static_controller<false>::problem_dimension);
    const ctrlpp_test::scripted_answer finite{ctrlpp::solve_status::optimal};

    SECTION("a primal one entry short is refused and the stored primal is kept")
    {
        auto record = scripted({finite, ctrlpp_test::scripted_answer{ctrlpp::solve_status::optimal, ctrlpp_test::corrupted_field::none, ctrlpp_test::report_values::nan, -1}});
        auto controller = make_scripted_static<false>(record, nonlinear_config(shape_horizon));
        REQUIRE(controller.solve(Eigen::Vector2d{1.0, 0.0}).has_value());

        check_refused(controller);
        CHECK(controller.last_solution().size() == posed);
        CHECK(bitwise_equal(controller.last_solution(), record->primals[0]));
    }

    SECTION("a primal one entry long is consumed up to the posed dimension")
    {
        auto record = scripted({ctrlpp_test::scripted_answer{ctrlpp::solve_status::optimal, ctrlpp_test::corrupted_field::none, ctrlpp_test::report_values::nan, 1}});
        auto controller = make_scripted_static<false>(record, nonlinear_config(shape_horizon));
        REQUIRE(controller.solve(Eigen::Vector2d{1.0, 0.0}).has_value());

        CHECK(controller.last_solution().size() == posed);
        CHECK(bitwise_equal(controller.last_solution(), record->primals[0].head(posed)));
    }
}

namespace
{

using construction_outcome = std::optional<ctrlpp::controller_construction_error>;

const construction_outcome accepted{};

constexpr double infinity = std::numeric_limits<double>::infinity();
constexpr double not_a_number = std::numeric_limits<double>::quiet_NaN();
constexpr double denorm = std::numeric_limits<double>::denorm_min();

using constrained_nonlinear_config = ctrlpp::nmpc_config<double, NX, NU, 1, 1>;
using constrained_nonlinear_controller = ctrlpp::nmpc_dynamic<double, NX, NU, ctrlpp_test::stub_nlp_solver<double>, decltype(double_integrator), 1, 1>;
using static_nonlinear_controller = ctrlpp::nmpc_static<double, NX, NU, static_cast<std::size_t>(shape_horizon), ctrlpp_test::stub_nlp_solver<double>, decltype(double_integrator)>;

using linear_construction_outcome = std::optional<linear_controller::construction_error>;

auto linear_refusal(construction_outcome outcome) -> linear_construction_outcome
{
    if(!outcome)
        return std::nullopt;
    return linear_controller::construction_error{*outcome};
}

auto riccati_refusal(ctrlpp::dare_error cause) -> linear_construction_outcome
{
    return linear_controller::construction_error{cause};
}

template <typename Result>
auto agreed_refusal(const Result& two_argument, const Result& solver_taking) -> std::optional<typename Result::error_type>
{
    REQUIRE(two_argument.has_value() == solver_taking.has_value());
    if(two_argument.has_value())
        return std::nullopt;
    REQUIRE(two_argument.error() == solver_taking.error());
    return two_argument.error();
}

template <typename Result>
auto agreed_outcome(const Result& two_argument, const Result& solver_taking) -> construction_outcome
{
    auto const refusal = agreed_refusal(two_argument, solver_taking);
    if(!refusal)
        return std::nullopt;
    auto const* configuration = std::get_if<ctrlpp::controller_construction_error>(&*refusal);
    REQUIRE(configuration != nullptr);
    return *configuration;
}

template <typename Configure>
auto linear_outcome_on(const ctrlpp::discrete_state_space<double, NX, NU, NX>& system, Configure configure) -> linear_construction_outcome
{
    auto cfg = linear_config(shape_horizon);
    configure(cfg);
    return agreed_refusal(linear_controller::create(system, cfg), linear_controller::create(system, cfg, ctrlpp_test::stub_qp_solver<double>{}));
}

template <typename Configure>
auto linear_outcome(Configure configure) -> linear_construction_outcome
{
    return linear_outcome_on(make_system(), configure);
}

template <typename Controller = nonlinear_controller, typename Config = ctrlpp::nmpc_config<double, NX, NU>, typename Configure>
auto nonlinear_outcome(Configure configure) -> construction_outcome
{
    Config cfg;
    cfg.horizon = shape_horizon;
    configure(cfg);
    return agreed_outcome(Controller::create(double_integrator, cfg), Controller::create(double_integrator, cfg, ctrlpp_test::stub_nlp_solver<double>{}));
}

template <typename Configure>
void require_both_predictive(Configure configure, construction_outcome expected)
{
    CHECK(linear_outcome(configure) == linear_refusal(expected));
    CHECK(nonlinear_outcome(configure) == expected);
}

struct input_bounds
{
    static constexpr auto refusal = ctrlpp::controller_construction_error::invalid_input_bounds;

    template <typename Config>
    static void set(Config& cfg, std::optional<double> lower, std::optional<double> upper)
    {
        if(lower)
            cfg.u_min = Eigen::Matrix<double, 1, 1>{*lower};
        if(upper)
            cfg.u_max = Eigen::Matrix<double, 1, 1>{*upper};
    }
};

// The first component is a valid box, so a refusal can only come from the
// second, which is the one under test.
struct state_bounds
{
    static constexpr auto refusal = ctrlpp::controller_construction_error::invalid_state_bounds;

    template <typename Config>
    static void set(Config& cfg, std::optional<double> lower, std::optional<double> upper)
    {
        if(lower)
            cfg.x_min = Eigen::Vector2d{-1.0, *lower};
        if(upper)
            cfg.x_max = Eigen::Vector2d{1.0, *upper};
    }
};

template <typename Bounds>
void require_bound_rules()
{
    auto check = [](std::optional<double> lower, std::optional<double> upper, construction_outcome expected)
    {
        INFO("lower " << lower.value_or(0.0) << (lower ? "" : " (unset)") << ", upper " << upper.value_or(0.0) << (upper ? "" : " (unset)"));
        require_both_predictive([=](auto& cfg) { Bounds::set(cfg, lower, upper); }, expected);
    };

    check(not_a_number, std::nullopt, Bounds::refusal);
    check(std::nullopt, not_a_number, Bounds::refusal);
    check(not_a_number, 1.0, Bounds::refusal);
    check(-1.0, not_a_number, Bounds::refusal);
    check(infinity, std::nullopt, Bounds::refusal);
    check(std::nullopt, -infinity, Bounds::refusal);
    check(1.0, -1.0, Bounds::refusal);
    check(1.0, std::nextafter(1.0, -infinity), Bounds::refusal);

    check(1.0, 1.0, accepted);
    check(-infinity, std::nullopt, accepted);
    check(std::nullopt, infinity, accepted);
    check(-infinity, infinity, accepted);
    check(-infinity, -1.0, accepted);
    check(1.0, infinity, accepted);
}

auto rate_bound(double value)
{
    return [value](auto& cfg) { cfg.du_max = Eigen::Matrix<double, 1, 1>{value}; };
}

// The symmetry bound and the definiteness floor restate the library's
// rounding_bound(N, scale, 1) = N * eps * scale at a scale where its
// denorm_min term is inactive, so a change to either derivation is caught.
auto rounding_bound_2x2(double scale) -> double
{
    return 2.0 * std::numeric_limits<double>::epsilon() * scale;
}

auto weight_2x2(double a, double b, double c, double d) -> Eigen::Matrix2d
{
    Eigen::Matrix2d weight;
    weight << a, b, c, d;
    return weight;
}

auto state_weight(const Eigen::Matrix2d& weight)
{
    return [weight](auto& cfg) { cfg.Q = weight; };
}

auto terminal_weight(const Eigen::Matrix2d& weight)
{
    return [weight](auto& cfg) { cfg.Qf = weight; };
}

auto input_weight(double weight)
{
    return [weight](auto& cfg) { cfg.R = Eigen::Matrix<double, 1, 1>{weight}; };
}

template <typename Weight>
void require_square_weight_rules(Weight weight)
{
    const double past_symmetry = std::nextafter(rounding_bound_2x2(2.0), infinity);
    const double past_floor = std::nextafter(rounding_bound_2x2(1.0), infinity);

    require_both_predictive(weight(weight_2x2(not_a_number, 0.0, 0.0, 1.0)), ctrlpp::controller_construction_error::non_finite_cost_weight);
    require_both_predictive(weight(weight_2x2(1.0, infinity, infinity, 1.0)), ctrlpp::controller_construction_error::non_finite_cost_weight);
    require_both_predictive(weight(weight_2x2(1.0, 0.5, 0.4, 1.0)), ctrlpp::controller_construction_error::asymmetric_cost_weight);
    require_both_predictive(weight(weight_2x2(2.0, 0.0, past_symmetry, 1.0)), ctrlpp::controller_construction_error::asymmetric_cost_weight);
    require_both_predictive(weight(weight_2x2(1.0, 0.0, 0.0, -1.0)), ctrlpp::controller_construction_error::indefinite_cost_weight);
    require_both_predictive(weight(weight_2x2(1.0, 0.0, 0.0, -past_floor)), ctrlpp::controller_construction_error::indefinite_cost_weight);

    require_both_predictive(weight(weight_2x2(2.0, 0.0, rounding_bound_2x2(2.0), 1.0)), accepted);
    require_both_predictive(weight(weight_2x2(2.0, 0.5, 0.5 + (std::nextafter(2.0, infinity) - 2.0), 1.0)), accepted);
    require_both_predictive(weight(weight_2x2(1.0, 0.0, 0.0, -rounding_bound_2x2(1.0))), accepted);
    require_both_predictive(weight(weight_2x2(1.0, 1.0, 1.0, 1.0)), accepted);
}

auto path_constrained(double penalty, bool soft)
{
    return [=](constrained_nonlinear_config& cfg)
    {
        cfg.soft_constraints = soft;
        cfg.path_constraint = [](const Eigen::Vector2d&, const Eigen::Matrix<double, 1, 1>&) -> Eigen::Matrix<double, 1, 1> { return Eigen::Matrix<double, 1, 1>::Zero(); };
        cfg.path_penalty << penalty;
    };
}

auto terminal_constrained(double penalty, bool soft)
{
    return [=](constrained_nonlinear_config& cfg)
    {
        cfg.soft_constraints = soft;
        cfg.terminal_constraint = [](const Eigen::Vector2d&) -> Eigen::Matrix<double, 1, 1> { return Eigen::Matrix<double, 1, 1>::Zero(); };
        cfg.terminal_penalty << penalty;
    };
}

auto unconstrained_penalties(double penalty)
{
    return [=](constrained_nonlinear_config& cfg)
    {
        cfg.path_penalty << penalty;
        cfg.terminal_penalty << penalty;
    };
}

auto constrained_outcome(auto configure) -> construction_outcome
{
    return nonlinear_outcome<constrained_nonlinear_controller, constrained_nonlinear_config>(configure);
}

auto soft_state_box(double penalty, bool hard, bool per_state)
{
    return [=](ctrlpp::mpc_config<double, NX, NU>& cfg)
    {
        cfg.x_min = Eigen::Vector2d{-1.0, -1.0};
        cfg.x_max = Eigen::Vector2d{1.0, 1.0};
        cfg.hard_state_constraints = hard;
        if(per_state)
            cfg.soft_state_penalty = Eigen::Vector2d{1.0, penalty};
        else
            cfg.soft_penalty = penalty;
    };
}

auto unbounded_state_penalty(double penalty, bool per_state)
{
    return [=](ctrlpp::mpc_config<double, NX, NU>& cfg)
    {
        if(per_state)
            cfg.soft_state_penalty = Eigen::Vector2d{1.0, penalty};
        else
            cfg.soft_penalty = penalty;
    };
}

const std::array<double, 6> invalid_penalties{0.0, -0.0, -denorm, -1.0, not_a_number, infinity};
const std::array<double, 2> valid_penalties{denorm, 1.0};

}

TEST_CASE("Predictive factories refuse inconsistent or non-finite input bounds", "[mpc][nmpc][configuration][hardening]")
{
    require_bound_rules<input_bounds>();
}

TEST_CASE("Predictive factories refuse inconsistent or non-finite state bounds", "[mpc][nmpc][configuration][hardening]")
{
    require_bound_rules<state_bounds>();
}

TEST_CASE("Predictive factories refuse a NaN or negative rate bound", "[mpc][nmpc][configuration][hardening]")
{
    require_both_predictive(rate_bound(not_a_number), ctrlpp::controller_construction_error::invalid_rate_bound);
    require_both_predictive(rate_bound(-denorm), ctrlpp::controller_construction_error::invalid_rate_bound);
    require_both_predictive(rate_bound(-infinity), ctrlpp::controller_construction_error::invalid_rate_bound);

    require_both_predictive(rate_bound(0.0), accepted);
    require_both_predictive(rate_bound(-0.0), accepted);
    require_both_predictive(rate_bound(denorm), accepted);
    require_both_predictive(rate_bound(infinity), accepted);
}

TEST_CASE("Predictive factories refuse a non-finite, structurally asymmetric or indefinite state weight", "[mpc][nmpc][configuration][hardening]")
{
    require_square_weight_rules([](const Eigen::Matrix2d& weight) { return state_weight(weight); });

    // A zero state weight poses a convex program, but no stabilizing Riccati solution exists for it on
    // this plant, so only the linear MPC with Qf unset refuses it, with the solve's own cause.
    auto const zero_state_weight = state_weight(Eigen::Matrix2d::Zero());
    CHECK(nonlinear_outcome(zero_state_weight) == accepted);
    CHECK(nonlinear_outcome<static_nonlinear_controller>(zero_state_weight) == accepted);
    CHECK(linear_outcome([&](auto& cfg) {
              zero_state_weight(cfg);
              cfg.Qf = Eigen::Matrix2d::Identity();
          })
          == linear_refusal(accepted));
    CHECK(linear_outcome(zero_state_weight) == riccati_refusal(ctrlpp::dare_error::non_stabilizable));
}

TEST_CASE("Predictive factories refuse a non-finite, structurally asymmetric or indefinite terminal weight when one is set", "[mpc][nmpc][configuration][hardening]")
{
    require_square_weight_rules([](const Eigen::Matrix2d& weight) { return terminal_weight(weight); });
    require_both_predictive(terminal_weight(Eigen::Matrix2d::Zero()), accepted);
    require_both_predictive([](auto& cfg) { cfg.Qf.reset(); }, accepted);
}

TEST_CASE("Predictive factories refuse a non-finite or negative input weight", "[mpc][nmpc][configuration][hardening]")
{
    require_both_predictive(input_weight(not_a_number), ctrlpp::controller_construction_error::non_finite_cost_weight);
    require_both_predictive(input_weight(-infinity), ctrlpp::controller_construction_error::non_finite_cost_weight);
    require_both_predictive(input_weight(-1.0), ctrlpp::controller_construction_error::indefinite_cost_weight);

    require_both_predictive(input_weight(1.0), accepted);
}

// The QP needs only a convex Hessian, so a singular R is refused only where the
// Riccati solve that supplies a missing terminal weight has to invert it.
TEST_CASE("A singular input weight is refused only where the Riccati solve inverts it", "[mpc][nmpc][configuration][hardening]")
{
    CHECK(linear_outcome(input_weight(0.0)) == linear_refusal(ctrlpp::controller_construction_error::singular_input_weight));
    CHECK(linear_outcome([](auto& cfg) {
              input_weight(0.0)(cfg);
              cfg.Qf = 10.0 * Eigen::Matrix2d::Identity();
          })
          == linear_refusal(accepted));

    CHECK(nonlinear_outcome(input_weight(0.0)) == accepted);
    CHECK(nonlinear_outcome<static_nonlinear_controller>(input_weight(0.0)) == accepted);
}

namespace
{

template <typename Corrupt>
auto corrupted_system(Corrupt corrupt) -> ctrlpp::discrete_state_space<double, NX, NU, NX>
{
    auto system = make_system();
    corrupt(system);
    return system;
}

auto scaled_output_map(double scale) -> ctrlpp::discrete_state_space<double, NX, NU, NX>
{
    return corrupted_system([scale](auto& system) { system.C *= scale; });
}

auto nan_dynamics() -> ctrlpp::discrete_state_space<double, NX, NU, NX>
{
    return corrupted_system([](auto& system) { system.A(0, 1) = not_a_number; });
}

void leave_configuration(ctrlpp::mpc_config<double, NX, NU>&) {}

void require_plant_outcome(const ctrlpp::discrete_state_space<double, NX, NU, NX>& system, const linear_construction_outcome& expected)
{
    CHECK(linear_outcome_on(system, leave_configuration) == expected);
    CHECK(linear_outcome_on(system, terminal_weight(Eigen::Matrix2d::Identity())) == expected);
}

}

TEST_CASE("Linear MPC refuses a non-finite plant with or without a terminal weight", "[mpc][configuration][hardening]")
{
    auto const refused = linear_refusal(ctrlpp::controller_construction_error::non_finite_system);
    require_plant_outcome(nan_dynamics(), refused);
    require_plant_outcome(corrupted_system([](auto& system) { system.A(1, 1) = -infinity; }), refused);
    require_plant_outcome(corrupted_system([](auto& system) { system.B(1) = infinity; }), refused);
    require_plant_outcome(corrupted_system([](auto& system) { system.C(1, 0) = not_a_number; }), refused);

    // The controller never reads D.
    require_plant_outcome(corrupted_system([](auto& system) { system.D(0) = not_a_number; }), linear_refusal(accepted));
}

TEST_CASE("Linear MPC refuses a state-space weight that overflows from finite operands", "[mpc][configuration][hardening]")
{
    // C = 1e200 I maps Q = I to 1e400 I, which overflows. C = 1e150 I keeps C^T Q C = 1e300 I representable,
    // and then only Qf = 1e10 I overflows.
    auto const refused = linear_refusal(ctrlpp::controller_construction_error::non_finite_state_space_weight);
    require_plant_outcome(scaled_output_map(1e200), refused);
    CHECK(linear_outcome_on(scaled_output_map(1e150), terminal_weight(1e10 * Eigen::Matrix2d::Identity())) == refused);
    CHECK(linear_outcome_on(scaled_output_map(1e150), terminal_weight(Eigen::Matrix2d::Identity())) == linear_refusal(accepted));
}

// Admission accepts q22 = -2 eps, exactly its floor. The second mode is stable at 0.8125 and uncontrollable, so the exact
// stabilizing solution has p22 = q22 / (1 - 0.8125^2) = -512 eps / 87, beyond the Riccati check's band of about 4.02 eps.
TEST_CASE("Linear MPC with Qf unset refuses a state weight admitted at its floor once the Riccati solution amplifies it", "[mpc][nmpc][configuration][hardening]")
{
    auto plant = make_system();
    plant.A = Eigen::Vector2d{0.1, 0.8125}.asDiagonal();
    plant.B = Eigen::Vector2d{1.0, 0.0};
    const Eigen::Matrix2d at_floor = weight_2x2(1.0, 0.0, 0.0, -rounding_bound_2x2(1.0));
    auto const weight = state_weight(at_floor);

    CHECK(nonlinear_outcome(weight) == accepted);
    CHECK(nonlinear_outcome<static_nonlinear_controller>(weight) == accepted);
    CHECK(linear_outcome_on(plant, [&](auto& cfg) {
              weight(cfg);
              cfg.Qf = Eigen::Matrix2d::Identity();
          })
          == linear_refusal(accepted));

    auto const riccati = ctrlpp::dare<double, NX, NU>(plant.A, plant.B, at_floor, Eigen::Matrix<double, 1, 1>::Identity());
    REQUIRE_FALSE(riccati.has_value());
    CHECK(riccati.error() == ctrlpp::dare_error::non_psd_solution);
    CHECK(linear_outcome_on(plant, weight) == riccati_refusal(ctrlpp::dare_error::non_psd_solution));
}

TEST_CASE("Runtime-horizon nonlinear MPC poses one symmetric weight when the configured one drifts by rounding", "[nmpc][configuration][hardening]")
{
    const double drift = std::nextafter(2.0, infinity) - 2.0;
    auto cfg = nonlinear_config(shape_horizon);
    cfg.Q = weight_2x2(2.0, 0.5, 0.5 + drift, 1.0);
    REQUIRE(cfg.Q(0, 1) != cfg.Q(1, 0));
    const Eigen::Matrix2d symmetric = weight_2x2(2.0, 0.5 + drift / 2.0, 0.5 + drift / 2.0, 1.0);

    auto controller = nonlinear_controller::create(double_integrator, cfg, ctrlpp_test::stub_nlp_solver<double>{});
    REQUIRE(controller.has_value());
    const auto& problem = controller->problem();

    std::vector<double> z(static_cast<std::size_t>(problem.n_vars), 0.0);
    std::vector<double> gradient(z.size(), 0.0);
    Eigen::Matrix2d posed;
    for(int column = 0; column < 2; ++column)
    {
        std::fill(z.begin(), z.end(), 0.0);
        z[static_cast<std::size_t>(column)] = 1.0;
        problem.gradient_cb(z, gradient);
        posed.col(column) = Eigen::Vector2d{gradient[0], gradient[1]};
    }
    CHECK((posed.array() == symmetric.array()).all());

    std::fill(z.begin(), z.end(), 0.0);
    z[0] = z[1] = 1.0;
    const Eigen::Vector2d dx{1.0, 1.0};
    CHECK(problem.cost_cb(z) == 0.5 * dx.dot(symmetric * dx));
}

TEST_CASE("Runtime-horizon nonlinear MPC refuses an invalid penalty only while its soft constraint is active", "[nmpc][configuration][hardening]")
{
    for(const double penalty : invalid_penalties)
    {
        INFO("penalty " << penalty);
        CHECK(constrained_outcome(path_constrained(penalty, true)) == ctrlpp::controller_construction_error::invalid_constraint_penalty);
        CHECK(constrained_outcome(terminal_constrained(penalty, true)) == ctrlpp::controller_construction_error::invalid_constraint_penalty);

        CHECK(constrained_outcome(path_constrained(penalty, false)) == accepted);
        CHECK(constrained_outcome(terminal_constrained(penalty, false)) == accepted);
        CHECK(constrained_outcome(unconstrained_penalties(penalty)) == accepted);
    }
    for(const double penalty : valid_penalties)
    {
        INFO("penalty " << penalty);
        CHECK(constrained_outcome(path_constrained(penalty, true)) == accepted);
        CHECK(constrained_outcome(terminal_constrained(penalty, true)) == accepted);
    }
}

TEST_CASE("Linear MPC refuses an invalid slack penalty only while its soft state constraints are active", "[mpc][configuration][hardening]")
{
    for(const bool per_state : {false, true})
    {
        for(const double penalty : invalid_penalties)
        {
            INFO("penalty " << penalty << (per_state ? " per state" : " uniform"));
            CHECK(linear_outcome(soft_state_box(penalty, false, per_state)) == linear_refusal(ctrlpp::controller_construction_error::invalid_constraint_penalty));

            CHECK(linear_outcome(soft_state_box(penalty, true, per_state)) == linear_refusal(accepted));
            CHECK(linear_outcome(unbounded_state_penalty(penalty, per_state)) == linear_refusal(accepted));
        }
        for(const double penalty : valid_penalties)
        {
            INFO("penalty " << penalty << (per_state ? " per state" : " uniform"));
            CHECK(linear_outcome(soft_state_box(penalty, false, per_state)) == linear_refusal(accepted));
        }
    }

    // The uniform penalty weighs nothing while a per-state one is set, but it is
    // still a configuration field with a valid default.
    CHECK(linear_outcome([](ctrlpp::mpc_config<double, NX, NU>& cfg) {
              soft_state_box(1.0, false, true)(cfg);
              cfg.soft_penalty = not_a_number;
          })
          == linear_refusal(ctrlpp::controller_construction_error::invalid_constraint_penalty));
}

namespace
{

template <typename Configure>
auto refusing_setup_outcome_on(const ctrlpp::discrete_state_space<double, NX, NU, NX>& system, Configure configure) -> linear_construction_outcome
{
    auto cfg = linear_config(shape_horizon);
    configure(cfg);
    ctrlpp_test::stub_qp_solver<double> solver;
    solver.fail_setup = true;
    auto refused = linear_controller::create(system, cfg, std::move(solver));
    if(refused.has_value())
        return std::nullopt;
    return refused.error();
}

template <typename Configure>
void require_refused_before_setup(const ctrlpp::discrete_state_space<double, NX, NU, NX>& system, Configure configure, const linear_construction_outcome& expected)
{
    CHECK(linear_outcome_on(system, configure) == expected);
    CHECK(refusing_setup_outcome_on(system, configure) == expected);
}

}

TEST_CASE("Linear MPC checks every configuration rule before setting its solver up", "[mpc][setup][configuration][hardening]")
{
    require_refused_before_setup(make_system(), [](auto& cfg) { input_bounds::set(cfg, 1.0, -1.0); }, linear_refusal(ctrlpp::controller_construction_error::invalid_input_bounds));
    require_refused_before_setup(nan_dynamics(), leave_configuration, linear_refusal(ctrlpp::controller_construction_error::non_finite_system));
    require_refused_before_setup(make_system(), input_weight(0.0), linear_refusal(ctrlpp::controller_construction_error::singular_input_weight));

    auto const plant = make_system();
    auto const riccati = ctrlpp::dare<double, NX, NU>(plant.A, plant.B, Eigen::Matrix2d::Zero(), Eigen::Matrix<double, 1, 1>::Identity());
    REQUIRE_FALSE(riccati.has_value());
    CHECK(riccati.error() == ctrlpp::dare_error::non_stabilizable);
    require_refused_before_setup(plant, state_weight(Eigen::Matrix2d::Zero()), riccati_refusal(riccati.error()));
}

namespace
{

using halfplanes = Eigen::Matrix<double, Eigen::Dynamic, 2>;

auto ellipsoid(const Eigen::Matrix2d& P, double alpha)
{
    return [=](ctrlpp::mpc_config<double, NX, NU>& cfg) { cfg.terminal_constraint_set = ctrlpp::ellipsoidal_set<double, NX>{P, alpha}; };
}

auto polytope(const halfplanes& H, const Eigen::VectorXd& h)
{
    return [=](ctrlpp::mpc_config<double, NX, NU>& cfg) { cfg.terminal_constraint_set = ctrlpp::polytopic_set<double, NX>{H, h}; };
}

auto box_rows() -> halfplanes
{
    halfplanes H(4, 2);
    H << 1.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, -1.0;
    return H;
}

auto box_levels() -> Eigen::VectorXd
{
    return Eigen::VectorXd::Constant(4, 0.08);
}

auto box_with_level(Eigen::Index row, double level)
{
    Eigen::VectorXd h = box_levels();
    h(row) = level;
    return polytope(box_rows(), h);
}

auto nan_level_box()
{
    Eigen::VectorXd h = box_levels();
    h(0) = h(2) = not_a_number;
    return polytope(box_rows(), h);
}

auto refused_terminal(ctrlpp::controller_construction_error cause) -> linear_construction_outcome
{
    return linear_refusal(cause);
}

template <typename Configure>
auto posed_problem(Configure configure) -> ctrlpp::qp_problem<double>
{
    REQUIRE(linear_outcome(configure) == linear_refusal(accepted));
    auto cfg = linear_config(shape_horizon);
    configure(cfg);
    ctrlpp::qp_problem<double> problem;
    ctrlpp_test::stub_qp_solver<double> solver;
    solver.posed = std::ref(problem);
    REQUIRE(linear_controller::create(make_system(), cfg, std::move(solver)).has_value());
    return problem;
}

template <typename Configure>
void require_posed_without_nan(Configure configure)
{
    auto const problem = posed_problem(configure);
    CHECK_FALSE(problem.l.array().isNaN().any());
    CHECK_FALSE(problem.u.array().isNaN().any());
    CHECK((problem.l.array() <= problem.u.array()).all());
}

auto dense_entries(const Eigen::SparseMatrix<double>& sparse) -> Eigen::VectorXd
{
    const Eigen::MatrixXd dense{sparse};
    return Eigen::Map<const Eigen::VectorXd>(dense.data(), dense.size());
}

auto random_unit_interval(std::mt19937_64& generator) -> double
{
    return 2.0 * (static_cast<double>(generator() >> 11) * 0x1.0p-53) - 1.0;
}

}

TEST_CASE("Linear MPC refuses a terminal set it cannot pose", "[mpc][configuration][terminal_set][hardening]")
{
    using enum ctrlpp::controller_construction_error;
    const double past_symmetry = std::nextafter(rounding_bound_2x2(2.0), infinity);
    const double past_floor = std::nextafter(rounding_bound_2x2(1.0), infinity);
    const Eigen::Matrix2d identity = Eigen::Matrix2d::Identity();

    CHECK(linear_outcome(ellipsoid(weight_2x2(not_a_number, 0.0, 0.0, 1.0), 1.0)) == refused_terminal(non_finite_terminal_set));
    CHECK(linear_outcome(ellipsoid(weight_2x2(1.0, infinity, infinity, 1.0), 1.0)) == refused_terminal(non_finite_terminal_set));
    CHECK(linear_outcome(ellipsoid(identity, not_a_number)) == refused_terminal(non_finite_terminal_set));
    CHECK(linear_outcome(ellipsoid(identity, infinity)) == refused_terminal(non_finite_terminal_set));
    CHECK(linear_outcome(ellipsoid(weight_2x2(1.0, 0.5, 0.4, 1.0), 1.0)) == refused_terminal(asymmetric_terminal_set));
    CHECK(linear_outcome(ellipsoid(weight_2x2(2.0, 0.0, past_symmetry, 1.0), 1.0)) == refused_terminal(asymmetric_terminal_set));
    CHECK(linear_outcome(ellipsoid(weight_2x2(1.0, 0.0, 0.0, -1.0), 1.0)) == refused_terminal(indefinite_terminal_set));
    CHECK(linear_outcome(ellipsoid(weight_2x2(1.0, 0.0, 0.0, -past_floor), 1.0)) == refused_terminal(indefinite_terminal_set));
    CHECK(linear_outcome(ellipsoid(identity, -1.0)) == refused_terminal(negative_terminal_set_level));
    CHECK(linear_outcome(ellipsoid(identity, -denorm)) == refused_terminal(negative_terminal_set_level));

    CHECK(linear_outcome(polytope(box_rows(), Eigen::VectorXd::Constant(3, 0.08))) == refused_terminal(mismatched_terminal_set));
    halfplanes nan_rows = box_rows();
    nan_rows(1, 0) = not_a_number;
    CHECK(linear_outcome(polytope(nan_rows, box_levels())) == refused_terminal(non_finite_terminal_set));
    halfplanes infinite_rows = box_rows();
    infinite_rows(2, 1) = infinity;
    CHECK(linear_outcome(polytope(infinite_rows, box_levels())) == refused_terminal(non_finite_terminal_set));
    CHECK(linear_outcome(nan_level_box()) == refused_terminal(non_finite_terminal_set));
}

TEST_CASE("Linear MPC poses every terminal set it admits with no NaN bound", "[mpc][configuration][terminal_set][hardening]")
{
    const Eigen::Matrix2d singular = Eigen::Vector2d{1.0, 0.0}.asDiagonal();
    const Eigen::Matrix2d at_floor = weight_2x2(1.0, 0.0, 0.0, -rounding_bound_2x2(1.0));
    for(const double alpha : {1.0, 0.0})
    {
        INFO("alpha " << alpha);
        require_posed_without_nan(ellipsoid(singular, alpha));
        require_posed_without_nan(ellipsoid(at_floor, alpha));
    }
    require_posed_without_nan(ellipsoid(Eigen::Matrix2d::Identity(), 0.0));
    require_posed_without_nan(ellipsoid(Eigen::Matrix2d::Zero(), 0.0));

    std::mt19937_64 generator{1};
    bool negative_eigenvalue_seen = false;
    for(int draw = 0; draw < 64; ++draw)
    {
        const Eigen::Vector2d v{random_unit_interval(generator), random_unit_interval(generator)};
        const Eigen::Matrix2d rank_one = v * v.transpose();
        negative_eigenvalue_seen = negative_eigenvalue_seen || Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d>{rank_one}.eigenvalues()(0) < 0.0;
        INFO("draw " << draw << ", v = (" << v(0) << ", " << v(1) << ")");
        require_posed_without_nan(ellipsoid(rank_one, 0.0));
        require_posed_without_nan(ellipsoid(rank_one, 1.0));
    }
    REQUIRE(negative_eigenvalue_seen);

    require_posed_without_nan(box_with_level(0, infinity));

    // Emptiness is not tested at construction: the row poses u = -inf, which the backend refuses at setup.
    auto const empty_row = posed_problem(box_with_level(0, -infinity));
    const Eigen::Index terminal_row = empty_row.u.size() - 4;
    CHECK_FALSE(empty_row.l.array().isNaN().any());
    CHECK_FALSE(empty_row.u.array().isNaN().any());
    CHECK(empty_row.u(terminal_row) == -infinity);
    CHECK((empty_row.l.head(terminal_row).array() <= empty_row.u.head(terminal_row).array()).all());
    CHECK((empty_row.l.tail(3).array() <= empty_row.u.tail(3).array()).all());
}

TEST_CASE("Linear MPC poses one symmetric ellipsoid weight when the configured one drifts by rounding", "[mpc][configuration][terminal_set][hardening]")
{
    const double drift = std::nextafter(2.0, infinity) - 2.0;
    const Eigen::Matrix2d drifted = weight_2x2(2.0, 0.5, 0.5 + drift, 1.0);
    REQUIRE(drifted(0, 1) != drifted(1, 0));
    const Eigen::Matrix2d symmetric = weight_2x2(2.0, 0.5 + drift / 2.0, 0.5 + drift / 2.0, 1.0);

    auto const posed = posed_problem(ellipsoid(drifted, 1.0));

    auto const plant = make_system();
    auto encoded = [&](const Eigen::Matrix2d& P)
    {
        const std::optional<ctrlpp::terminal_set<double, NX>> set{ctrlpp::ellipsoidal_set<double, NX>{P, 1.0}};
        auto A = ctrlpp::detail::build_constraint_matrix<double, NX, NU>(shape_horizon, plant.A, plant.B, std::nullopt, std::nullopt, false, std::nullopt, std::nullopt, std::nullopt, set);
        auto [l, u] = ctrlpp::detail::build_bounds_vectors<double, NX, NU>(shape_horizon, Eigen::Vector2d::Zero(), std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, set);
        return std::tuple{dense_entries(A), l, u};
    };
    auto const [A_symmetric, l_symmetric, u_symmetric] = encoded(symmetric);
    auto const [A_raw, l_raw, u_raw] = encoded(drifted);
    REQUIRE_FALSE((bitwise_equal(A_raw, A_symmetric) && bitwise_equal(l_raw, l_symmetric) && bitwise_equal(u_raw, u_symmetric)));

    CHECK(bitwise_equal(dense_entries(posed.A), A_symmetric));
    CHECK(bitwise_equal(posed.l, l_symmetric));
    CHECK(bitwise_equal(posed.u, u_symmetric));
}

TEST_CASE("Predictive factories report the first violated rule in the documented order", "[mpc][nmpc][configuration][hardening]")
{
    auto both = [](auto first, auto second) {
        return [=](auto& cfg)
        {
            first(cfg);
            second(cfg);
        };
    };
    auto inverted_input = [](auto& cfg) { input_bounds::set(cfg, 1.0, -1.0); };
    auto inverted_state = [](auto& cfg) { state_bounds::set(cfg, 1.0, -1.0); };
    auto nan_state_weight = state_weight(weight_2x2(not_a_number, 0.0, 0.0, 1.0));
    auto indefinite_state_weight = state_weight(weight_2x2(1.0, 0.0, 0.0, -1.0));

    require_both_predictive(both(inverted_state, inverted_input), ctrlpp::controller_construction_error::invalid_input_bounds);
    require_both_predictive(both(rate_bound(-1.0), inverted_state), ctrlpp::controller_construction_error::invalid_state_bounds);
    require_both_predictive(both(nan_state_weight, rate_bound(-1.0)), ctrlpp::controller_construction_error::invalid_rate_bound);
    require_both_predictive(both(input_weight(not_a_number), indefinite_state_weight), ctrlpp::controller_construction_error::indefinite_cost_weight);
    require_both_predictive(both(terminal_weight(weight_2x2(1.0, 0.5, 0.4, 1.0)), input_weight(-1.0)), ctrlpp::controller_construction_error::indefinite_cost_weight);

    CHECK(constrained_outcome(both(path_constrained(0.0, true), nan_state_weight)) == ctrlpp::controller_construction_error::non_finite_cost_weight);
    CHECK(linear_outcome(both(soft_state_box(0.0, false, false), nan_state_weight)) == linear_refusal(ctrlpp::controller_construction_error::non_finite_cost_weight));
    CHECK(linear_outcome(both(input_weight(0.0), inverted_input)) == linear_refusal(ctrlpp::controller_construction_error::invalid_input_bounds));
    CHECK(linear_outcome(both(input_weight(0.0), soft_state_box(0.0, false, false))) == linear_refusal(ctrlpp::controller_construction_error::invalid_constraint_penalty));

    auto zero_state_weight = state_weight(Eigen::Matrix2d::Zero());
    CHECK(linear_outcome(both(zero_state_weight, inverted_input)) == linear_refusal(ctrlpp::controller_construction_error::invalid_input_bounds));
    CHECK(linear_outcome(both(zero_state_weight, input_weight(0.0))) == linear_refusal(ctrlpp::controller_construction_error::singular_input_weight));

    CHECK(linear_outcome_on(nan_dynamics(), inverted_input) == linear_refusal(ctrlpp::controller_construction_error::invalid_input_bounds));
    CHECK(linear_outcome_on(nan_dynamics(), input_weight(0.0)) == linear_refusal(ctrlpp::controller_construction_error::non_finite_system));
    CHECK(linear_outcome_on(nan_dynamics(), zero_state_weight) == linear_refusal(ctrlpp::controller_construction_error::non_finite_system));
    auto const nan_dynamics_overflowing_map = corrupted_system([](auto& system) {
        system.A(0, 1) = not_a_number;
        system.C *= 1e200;
    });
    CHECK(linear_outcome_on(nan_dynamics_overflowing_map, leave_configuration) == linear_refusal(ctrlpp::controller_construction_error::non_finite_system));
    CHECK(linear_outcome_on(scaled_output_map(1e200), input_weight(0.0)) == linear_refusal(ctrlpp::controller_construction_error::non_finite_state_space_weight));

    auto const nan_box = nan_level_box();
    auto const non_finite_set = linear_refusal(ctrlpp::controller_construction_error::non_finite_terminal_set);
    CHECK(linear_outcome(both(nan_box, inverted_input)) == linear_refusal(ctrlpp::controller_construction_error::invalid_input_bounds));
    CHECK(linear_outcome(both(nan_box, soft_state_box(0.0, false, false))) == linear_refusal(ctrlpp::controller_construction_error::invalid_constraint_penalty));
    CHECK(linear_outcome_on(nan_dynamics(), nan_box) == non_finite_set);
    CHECK(refusing_setup_outcome_on(make_system(), nan_box) == non_finite_set);
    halfplanes nan_rows = box_rows();
    nan_rows(0, 0) = not_a_number;
    CHECK(linear_outcome(polytope(nan_rows, Eigen::VectorXd::Constant(3, 0.08))) == linear_refusal(ctrlpp::controller_construction_error::mismatched_terminal_set));
    CHECK(linear_outcome(ellipsoid(weight_2x2(1.0, 0.0, 0.0, -1.0), not_a_number)) == linear_refusal(ctrlpp::controller_construction_error::indefinite_terminal_set));

    auto negative_horizon = [](auto& cfg) { cfg.horizon = -1; };
    require_both_predictive(both(inverted_input, negative_horizon), ctrlpp::controller_construction_error::non_positive_horizon);
}
