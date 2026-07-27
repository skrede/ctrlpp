// This anchor asserts that an unconstrained linear MHE, driven by the same
// input/measurement sequence as an in-process Kalman filter on an identical
// linear plant, converges to the same state estimate once its window has
// filled. Full-information estimation with a correctly formed arrival cost
// and process-noise weighting reduces to the Kalman recursion for a linear
// system, so this equivalence is the ground truth.
//
// The process noise enters the estimate through the soft dynamics residual
// weighted by the inverse process covariance (there is no hard dynamics
// equality to nullify it), and the arrival cost is anchored at the window head
// using the filter's predicted estimate from before the window, so it
// summarizes the pre-window data without double-counting the measurements
// already inside the window. The window-end marginal of the resulting
// least-squares problem is the Kalman posterior at the current time.
//
// The estimator solves its quadratic program to a tight tolerance so the
// residual against the reference is set by the cross-method comparison (a
// splitting QP solver versus the filter's recursive orthogonal solves) rather
// than by the solver stopping early. The tolerance below is derived from that
// QP convergence tolerance and the solution scale, not from machine epsilon.

#include "ctrlpp/mhe.h"
#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/estimation/kalman.h"
#include "ctrlpp/model/state_space.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <limits>
#include <random>
#include <cstddef>

using namespace ctrlpp;
using Catch::Matchers::WithinAbs;

namespace
{

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr std::size_t NY = 1;
constexpr std::size_t N = 5;
constexpr double dt = 0.1;

struct linear_dynamics
{
    auto operator()(const Vector<double, NX>& x, const Vector<double, NU>& u) const -> Vector<double, NX>
    {
        Vector<double, NX> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    }

    auto jacobian_x(const Vector<double, NX>& /*x*/, const Vector<double, NU>& /*u*/) const -> Matrix<double, NX, NX>
    {
        Matrix<double, NX, NX> F;
        F << 1.0, dt, 0.0, 1.0;
        return F;
    }

    auto jacobian_u(const Vector<double, NX>& /*x*/, const Vector<double, NU>& /*u*/) const -> Matrix<double, NX, NU>
    {
        Matrix<double, NX, NU> G;
        G << 0.5 * dt * dt, dt;
        return G;
    }
};

struct position_measurement
{
    auto operator()(const Vector<double, NX>& x) const -> Vector<double, NY>
    {
        return (Vector<double, NY>() << x(0)).finished();
    }

    auto jacobian(const Vector<double, NX>& /*x*/) const -> Matrix<double, NY, NX>
    {
        return (Matrix<double, NY, NX>() << 1.0, 0.0).finished();
    }
};

// Convergence tolerance the moving-horizon QP is solved to. The equivalence is
// limited by how accurately the QP is solved, not by floating-point round-off,
// so the estimator here drives OSQP far below its loose default: at the default
// tolerance the operator-splitting iterate stops several orders short and the
// gap to the Kalman recursion is solver truncation, not a formulation error.
constexpr double eps_qp = 1e-10;

// Adapter that instantiates the QP backend at the tight tolerance above (with
// polishing requested). The estimator constructs its solver internally, so the
// tolerance is bound to the backend type rather than injected at call sites.
// For this constraint-free anchor the accuracy comes from the tight splitting
// tolerance; polishing is inactive because there is no active constraint set.
struct tight_qp_solver
{
    using scalar_type = double;

    osqp_solver inner_{eps_qp, eps_qp, 40000, false, true, true};

    auto setup(const qp_problem<double>& problem) -> ctrlpp::expected<void, osqp_setup_error> { return inner_.setup(problem); }
    auto solve(const qp_update<double>& update) -> qp_result<double> { return inner_.solve(update); }
};

using MheType = mhe<double, NX, NU, NY, N, tight_qp_solver, linear_dynamics, position_measurement>;

} // namespace

TEST_CASE("unconstrained linear MHE matches Kalman filter at the window end", "[mhe][kalman][anchor]")
{
    mhe_config<double, NX, NU, NY, N> mhe_cfg;
    mhe_cfg.Q = Matrix<double, NX, NX>::Identity() * 0.01;
    mhe_cfg.R = Matrix<double, NY, NY>::Identity() * 0.1;
    mhe_cfg.P0 = Matrix<double, NX, NX>::Identity() * 10.0;

    MheType estimator(linear_dynamics{}, position_measurement{}, mhe_cfg);

    discrete_state_space<double, NX, NU, NY> sys;
    sys.A << 1.0, dt, 0.0, 1.0;
    sys.B << 0.5 * dt * dt, dt;
    sys.C << 1.0, 0.0;
    sys.D = Matrix<double, NY, NU>::Zero();

    kalman_config<double, NX, NU, NY> kf_cfg;
    kf_cfg.Q = mhe_cfg.Q;
    kf_cfg.R = mhe_cfg.R;
    kf_cfg.P0 = mhe_cfg.P0;

    kalman_filter<double, NX, NU, NY> reference(sys, kf_cfg);

    Vector<double, NX> x_true;
    x_true << 0.0, 1.0;
    Vector<double, NU> u = Vector<double, NU>::Zero();

    std::mt19937 gen(42);
    std::normal_distribution<double> noise(0.0, 0.1);

    for(std::size_t i = 0; i < 3 * N; ++i)
    {
        x_true = linear_dynamics{}(x_true, u);
        Vector<double, NY> z;
        z << x_true(0) + noise(gen);

        estimator.predict(u);
        estimator.update(z);
        reference.predict(u);
        reference.update(z);
    }

    // These are two different numerical methods: the estimator solves a
    // quadratic program by operator splitting to a residual of eps_qp, while
    // the reference runs the Kalman recursion through exact orthogonal solves.
    // Two such methods cannot agree tighter than the QP convergence tolerance
    // times the solution magnitude, amplified by the coupling of the window's
    // stage blocks through which the per-solve residual propagates. The bound
    // is therefore N (the horizon length, i.e. the number of coupled stage
    // blocks) times eps_qp times the solution scale. This is deliberately not a
    // machine-epsilon bound: the limiting factor is the solver tolerance, which
    // is many orders above floating-point round-off, so scaling by epsilon here
    // would be far too tight for a cross-method comparison.
    const double scale = 1.0 + reference.state().norm();
    const double tol = static_cast<double>(N) * eps_qp * scale;

    REQUIRE_THAT(estimator.state()(0), WithinAbs(reference.state()(0), tol));
    REQUIRE_THAT(estimator.state()(1), WithinAbs(reference.state()(1), tol));
}
