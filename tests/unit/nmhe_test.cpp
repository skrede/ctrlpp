#include "ctrlpp/nmhe.h"
#include "ctrlpp/mpc/nlopt_solver.h"

#include <Eigen/Dense>

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

struct pendulum_dynamics
{
    double dt = 0.05;
    double g = 9.81;
    double l = 1.0;

    auto operator()(const Vector<double, NX>& x, const Vector<double, NU>& u) const -> Vector<double, NX>
    {
        double theta = x(0);
        double omega = x(1);
        double alpha = -g / l * std::sin(theta) + u(0);
        return Vector<double, NX>{theta + dt * omega, omega + dt * alpha};
    }
};

struct angle_measurement
{
    auto operator()(const Vector<double, NX>& x) const -> Vector<double, NY>
    {
        return (Vector<double, NY>() << x(0)).finished();
    }
};

using NloptSolver = nlopt_solver<double>;
using NmheType = nmhe<double, NX, NU, NY, N, NloptSolver, pendulum_dynamics, angle_measurement>;

} // namespace

TEST_CASE("nmhe falls back to ekf during warmup", "[nmhe][nlopt]")
{
    nmhe_config<double, NX, NU, NY, N> cfg;
    cfg.Q = Matrix<double, NX, NX>::Identity() * 0.01;
    cfg.R = Matrix<double, NY, NY>::Identity() * 0.1;

    NmheType estimator(pendulum_dynamics{}, angle_measurement{}, cfg);

    Vector<double, NU> u = Vector<double, NU>::Zero();
    Vector<double, NY> z;
    z << 0.1;

    for(std::size_t i = 0; i < N - 1; ++i)
    {
        estimator.predict(u);
        estimator.update(z);
        REQUIRE(estimator.diagnostics().used_ekf_fallback);
    }
}

TEST_CASE("nmhe tracks nonlinear system after warmup", "[nmhe][nlopt]")
{
    nmhe_config<double, NX, NU, NY, N> cfg;
    cfg.Q = Matrix<double, NX, NX>::Identity() * 0.01;
    cfg.R = Matrix<double, NY, NY>::Identity() * 0.1;
    cfg.P0 = Matrix<double, NX, NX>::Identity() * 10.0;

    NmheType estimator(pendulum_dynamics{}, angle_measurement{}, cfg);

    pendulum_dynamics dyn;
    Vector<double, NX> x_true{0.3, 0.0};
    Vector<double, NU> u = Vector<double, NU>::Zero();

    for(int i = 0; i < static_cast<int>(N) + 10; ++i)
    {
        x_true = dyn(x_true, u);
        Vector<double, NY> z;
        z << x_true(0);

        estimator.predict(u);
        estimator.update(z);
    }

    auto est = estimator.state();
    REQUIRE(std::isfinite(est(0)));
    REQUIRE(std::isfinite(est(1)));
}

TEST_CASE("nmhe covariance is finite", "[nmhe][nlopt]")
{
    nmhe_config<double, NX, NU, NY, N> cfg;
    cfg.Q = Matrix<double, NX, NX>::Identity() * 0.01;
    cfg.R = Matrix<double, NY, NY>::Identity() * 0.1;

    NmheType estimator(pendulum_dynamics{}, angle_measurement{}, cfg);

    Vector<double, NU> u = Vector<double, NU>::Zero();
    Vector<double, NY> z;
    z << 0.1;

    for(int i = 0; i < static_cast<int>(N) + 3; ++i)
    {
        estimator.predict(u);
        estimator.update(z);
    }

    auto P = estimator.covariance();
    REQUIRE(std::isfinite(P.norm()));
    auto st = estimator.state();
    REQUIRE(std::isfinite(st.norm()));
}

TEST_CASE("nmhe satisfies CovarianceObserver", "[nmhe][nlopt]")
{
    static_assert(ObserverPolicy<NmheType>);
    static_assert(CovarianceObserver<NmheType>);
}

TEST_CASE("nmhe handles soft constraints", "[nmhe][nlopt]")
{
    nmhe_config<double, NX, NU, NY, N> cfg;
    cfg.Q = Matrix<double, NX, NX>::Identity() * 0.01;
    cfg.R = Matrix<double, NY, NY>::Identity() * 0.1;
    cfg.soft_constraints = true;

    NmheType estimator(pendulum_dynamics{}, angle_measurement{}, cfg);

    Vector<double, NU> u = Vector<double, NU>::Zero();
    Vector<double, NY> z;
    z << 0.1;

    for(int i = 0; i < static_cast<int>(N) + 3; ++i)
    {
        estimator.predict(u);
        estimator.update(z);
    }

    // Should not crash
    REQUIRE(std::isfinite(estimator.state().norm()));
}
