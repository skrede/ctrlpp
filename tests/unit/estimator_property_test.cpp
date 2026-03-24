#include "ctrlpp/estimation/ekf.h"

#include <catch2/catch_test_macros.hpp>
#include <rapidcheck.h>
#include <rapidcheck/catch.h>

#include <Eigen/Eigenvalues>

#include <cmath>
#include <cstddef>

using namespace ctrlpp;

namespace
{

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr std::size_t NY = 1;

struct linear_dynamics
{
    double dt = 0.1;

    auto operator()(const Vector<double, NX>& x, const Vector<double, NU>& u) const -> Vector<double, NX>
    {
        Vector<double, NX> xn;
        xn(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        xn(1) = x(1) + dt * u(0);
        return xn;
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

auto bounded_double(double lo, double hi) -> rc::Gen<double>
{
    return rc::gen::map(rc::gen::inRange(0, 1000000), [lo, hi](int x) { return lo + (hi - lo) * (static_cast<double>(x) / 1000000.0); });
}

using EkfType = ekf<double, NX, NU, NY, linear_dynamics, position_measurement>;

} // namespace

TEST_CASE("estimator property tests", "[estimator][property]")
{
    SECTION("ekf covariance stays symmetric over random inputs")
    {
        rc::prop("covariance symmetry", [](void)
                 {
            auto q_diag = *bounded_double(0.001, 1.0);
            auto r_val = *bounded_double(0.01, 10.0);

            ekf_config<double, NX, NU, NY> cfg;
            cfg.Q = Matrix<double, NX, NX>::Identity() * q_diag;
            cfg.R = Matrix<double, NY, NY>::Identity() * r_val;
            cfg.P0 = Matrix<double, NX, NX>::Identity();

            EkfType filter(linear_dynamics{}, position_measurement{}, cfg);

            for(int i = 0; i < 30; ++i)
            {
                auto u_val = *bounded_double(-5.0, 5.0);
                auto z_val = *bounded_double(-10.0, 10.0);
                filter.predict((Vector<double, NU>() << u_val).finished());
                filter.update((Vector<double, NY>() << z_val).finished());
            }

            auto P = filter.covariance();
            double asym = (P - P.transpose()).norm();
            RC_ASSERT(asym < 1e-10); });
    }

    SECTION("ekf covariance stays PSD over random inputs")
    {
        rc::prop("covariance PSD", [](void)
                 {
            auto q_diag = *bounded_double(0.001, 1.0);
            auto r_val = *bounded_double(0.01, 10.0);

            ekf_config<double, NX, NU, NY> cfg;
            cfg.Q = Matrix<double, NX, NX>::Identity() * q_diag;
            cfg.R = Matrix<double, NY, NY>::Identity() * r_val;
            cfg.P0 = Matrix<double, NX, NX>::Identity();

            EkfType filter(linear_dynamics{}, position_measurement{}, cfg);

            for(int i = 0; i < 30; ++i)
            {
                auto u_val = *bounded_double(-5.0, 5.0);
                auto z_val = *bounded_double(-10.0, 10.0);
                filter.predict((Vector<double, NU>() << u_val).finished());
                filter.update((Vector<double, NY>() << z_val).finished());
            }

            auto P = filter.covariance();
            Eigen::SelfAdjointEigenSolver<Matrix<double, NX, NX>> es(P);
            RC_ASSERT(es.eigenvalues().minCoeff() >= -1e-10); });
    }

    SECTION("ekf state converges toward measurements")
    {
        rc::prop("state convergence", [](void)
                 {
            auto true_pos = *bounded_double(-5.0, 5.0);

            ekf_config<double, NX, NU, NY> cfg;
            cfg.Q = Matrix<double, NX, NX>::Identity() * 0.01;
            cfg.R = Matrix<double, NY, NY>::Identity() * 0.1;
            cfg.P0 = Matrix<double, NX, NX>::Identity() * 10.0;

            EkfType filter(linear_dynamics{}, position_measurement{}, cfg);

            double initial_err = std::abs(filter.state()(0) - true_pos);

            for(int i = 0; i < 50; ++i)
            {
                filter.predict(Vector<double, NU>::Zero());
                filter.update((Vector<double, NY>() << true_pos).finished());
            }

            double final_err = std::abs(filter.state()(0) - true_pos);
            RC_ASSERT(final_err < initial_err + 0.1); });
    }
}
