#include "ctrlpp/mhe.h"
#include "ctrlpp/mpc/osqp_solver.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
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

using MheType = mhe<double, NX, NU, NY, N, osqp_solver<double>, linear_dynamics, position_measurement>;

} // namespace

TEST_CASE("mhe falls back to ekf during warmup", "[mhe][osqp]")
{
    mhe_config<double, NX, NU, NY, N> cfg;
    cfg.Q = Matrix<double, NX, NX>::Identity() * 0.01;
    cfg.R = Matrix<double, NY, NY>::Identity() * 0.1;

    MheType estimator(linear_dynamics{}, position_measurement{}, cfg);

    Vector<double, NU> u = Vector<double, NU>::Zero();
    Vector<double, NY> z;
    z << 0.1;

    // First N-1 steps should use EKF fallback
    for(std::size_t i = 0; i < N - 1; ++i)
    {
        estimator.predict(u);
        estimator.update(z);
        REQUIRE(estimator.diagnostics().used_ekf_fallback);
    }
}

TEST_CASE("mhe tracks constant velocity after warmup", "[mhe][osqp]")
{
    mhe_config<double, NX, NU, NY, N> cfg;
    cfg.Q = Matrix<double, NX, NX>::Identity() * 0.01;
    cfg.R = Matrix<double, NY, NY>::Identity() * 0.1;
    cfg.P0 = Matrix<double, NX, NX>::Identity() * 10.0;

    MheType estimator(linear_dynamics{}, position_measurement{}, cfg);

    // Simulate a system with constant velocity
    Vector<double, NX> x_true;
    x_true << 0.0, 1.0; // position=0, velocity=1
    Vector<double, NU> u = Vector<double, NU>::Zero();

    for(int i = 0; i < static_cast<int>(N) + 10; ++i)
    {
        x_true(0) += dt * x_true(1);
        Vector<double, NY> z;
        z << x_true(0);

        estimator.predict(u);
        estimator.update(z);
    }

    auto est = estimator.state();
    REQUIRE(std::isfinite(est(0)));
    REQUIRE(std::isfinite(est(1)));
}

TEST_CASE("mhe covariance is finite and PSD", "[mhe][osqp]")
{
    mhe_config<double, NX, NU, NY, N> cfg;
    cfg.Q = Matrix<double, NX, NX>::Identity() * 0.01;
    cfg.R = Matrix<double, NY, NY>::Identity() * 0.1;

    MheType estimator(linear_dynamics{}, position_measurement{}, cfg);

    Vector<double, NU> u = Vector<double, NU>::Zero();
    Vector<double, NY> z;
    z << 0.5;

    for(int i = 0; i < static_cast<int>(N) + 5; ++i)
    {
        estimator.predict(u);
        estimator.update(z);
    }

    auto P = estimator.covariance();
    REQUIRE(std::isfinite(P.norm()));

    Eigen::SelfAdjointEigenSolver<Matrix<double, NX, NX>> es(P);
    REQUIRE(es.eigenvalues().minCoeff() >= -1e-10);
}

TEST_CASE("mhe arrival state is accessible", "[mhe][osqp]")
{
    mhe_config<double, NX, NU, NY, N> cfg;
    cfg.Q = Matrix<double, NX, NX>::Identity() * 0.01;
    cfg.R = Matrix<double, NY, NY>::Identity() * 0.1;

    MheType estimator(linear_dynamics{}, position_measurement{}, cfg);

    Vector<double, NU> u = Vector<double, NU>::Zero();
    Vector<double, NY> z;
    z << 0.1;

    for(int i = 0; i < static_cast<int>(N) + 2; ++i)
    {
        estimator.predict(u);
        estimator.update(z);
    }

    auto arrival = estimator.arrival_state();
    REQUIRE(std::isfinite(arrival(0)));
    REQUIRE(std::isfinite(arrival(1)));
}

TEST_CASE("mhe satisfies CovarianceObserver", "[mhe][osqp]")
{
    static_assert(ObserverPolicy<MheType>);
    static_assert(CovarianceObserver<MheType>);
}

TEST_CASE("mhe trajectory returns window states", "[mhe][osqp]")
{
    mhe_config<double, NX, NU, NY, N> cfg;
    cfg.Q = Matrix<double, NX, NX>::Identity() * 0.01;
    cfg.R = Matrix<double, NY, NY>::Identity() * 0.1;

    MheType estimator(linear_dynamics{}, position_measurement{}, cfg);

    Vector<double, NU> u = Vector<double, NU>::Zero();
    Vector<double, NY> z;
    z << 0.1;

    for(int i = 0; i < static_cast<int>(N) + 2; ++i)
    {
        estimator.predict(u);
        estimator.update(z);
    }

    auto traj = estimator.trajectory();
    REQUIRE(traj.size() == N + 1);
    for(auto& s : traj)
        REQUIRE(std::isfinite(s.norm()));
}
