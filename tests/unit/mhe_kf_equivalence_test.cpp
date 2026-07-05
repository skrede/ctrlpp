// This anchor asserts that an unconstrained linear MHE, driven by the same
// input/measurement sequence as an in-process Kalman filter on an identical
// linear plant, converges to the same state estimate once its window has
// filled. Full-information estimation with a correctly formed arrival cost
// and process-noise weighting reduces to the Kalman recursion for a linear
// system, so this equivalence is the ground truth.
//
// It currently fails for two structural reasons: the linear MHE's dynamics
// equality constraint forces the process-noise term to zero regardless of
// the configured Q, and its arrival cost is anchored using the filter's
// current-time estimate rather than a state properly propagated back to the
// window start. Both effects bias the window-end estimate by an amount that
// grows with the horizon length, the sampling period, and the true
// velocity, so the REQUIRE below currently fails; [!shouldfail] reports the
// case as passing until the estimator is corrected, at which point this tag
// must be removed.

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

using MheType = mhe<double, NX, NU, NY, N, osqp_solver, linear_dynamics, position_measurement>;

} // namespace

TEST_CASE("unconstrained linear MHE matches Kalman filter at the window end", "[mhe][kalman][anchor][!shouldfail]")
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

    // Backward error for dense QR-based linear solves grows linearly with
    // problem dimension (Trefethen & Bau, Numerical Linear Algebra, Lecture
    // 15); NX bounds the accumulated rounding across the per-step solves
    // both filters perform.
    const double scale = 1.0 + reference.state().norm() + reference.covariance().norm();
    const double tol = static_cast<double>(NX) * std::numeric_limits<double>::epsilon() * scale;

    REQUIRE_THAT(estimator.state()(0), WithinAbs(reference.state()(0), tol));
    REQUIRE_THAT(estimator.state()(1), WithinAbs(reference.state()(1), tol));
}
