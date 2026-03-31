#include "ctrlpp/estimation/particle_filter.h"
#include "ctrlpp/estimation/resampling/systematic_resampling.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstddef>
#include <random>

using namespace ctrlpp;

// ---------------------------------------------------------------------------
// Test dynamics: double integrator
// ---------------------------------------------------------------------------
struct pf_linear_dynamics
{
    double dt = 0.1;

    auto operator()(const Vector<double, 2>& x, const Vector<double, 1>& u) const -> Vector<double, 2>
    {
        Vector<double, 2> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    }
};

struct pf_position_measurement
{
    auto operator()(const Vector<double, 2>& x) const -> Vector<double, 1>
    {
        Vector<double, 1> z;
        z(0) = x(0);
        return z;
    }
};

// ---------------------------------------------------------------------------
// Core tracking tests
// ---------------------------------------------------------------------------

TEST_CASE("particle filter tracks linear system")
{
    pf_linear_dynamics dyn;
    pf_position_measurement meas;

    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.01;
    Matrix<double, 1, 1> R;
    R << 0.5;
    Vector<double, 2> x0 = Vector<double, 2>::Zero();
    Matrix<double, 2, 2> P0 = Matrix<double, 2, 2>::Identity() * 1.0;

    auto pf = make_particle_filter<500>(dyn, meas, pf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0}, std::mt19937_64{42});

    double true_pos = 0.0;
    double true_vel = 1.0;
    constexpr double dt = 0.1;

    for(int i = 0; i < 50; ++i)
    {
        true_pos += true_vel * dt;

        Vector<double, 1> u = Vector<double, 1>::Zero();
        pf.predict(u);

        Vector<double, 1> z;
        z << true_pos + 0.1 * std::sin(static_cast<double>(i));
        pf.update(z);
    }

    auto est = pf.state();
    CHECK_THAT(est(0), Catch::Matchers::WithinAbs(true_pos, 1.0));
    CHECK_THAT(est(1), Catch::Matchers::WithinAbs(true_vel, 1.0));
}

TEST_CASE("particle filter deterministic with fixed seed")
{
    pf_linear_dynamics dyn;
    pf_position_measurement meas;

    pf_config<double, 2, 1, 1> cfg{
        .Q = Matrix<double, 2, 2>::Identity() * 0.01,
        .R = ([]{ Matrix<double, 1, 1> R; R << 0.5; return R; })(),
        .x0 = Vector<double, 2>::Zero(),
        .P0 = Matrix<double, 2, 2>::Identity()
    };

    auto run = [&](std::uint64_t seed)
    {
        auto pf = make_particle_filter<100>(dyn, meas, cfg, std::mt19937_64{seed});

        Vector<double, 1> u = Vector<double, 1>::Zero();
        pf.predict(u);

        Vector<double, 1> z;
        z << 1.0;
        pf.update(z);

        return pf.state();
    };

    auto est1 = run(12345);
    auto est2 = run(12345);

    CHECK_THAT(est1(0), Catch::Matchers::WithinAbs(est2(0), 1e-14));
    CHECK_THAT(est1(1), Catch::Matchers::WithinAbs(est2(1), 1e-14));
}

TEST_CASE("particle filter weighted_mean and map_estimate both accessible")
{
    pf_linear_dynamics dyn;
    pf_position_measurement meas;

    auto pf = make_particle_filter<50>(dyn, meas, pf_config<double, 2, 1, 1>{}, std::mt19937_64{77});

    Vector<double, 1> u = Vector<double, 1>::Zero();
    pf.predict(u);

    Vector<double, 1> z;
    z << 1.0;
    pf.update(z);

    auto wm = pf.weighted_mean();
    auto me = pf.map_estimate();

    CHECK(std::isfinite(wm(0)));
    CHECK(std::isfinite(wm(1)));
    CHECK(std::isfinite(me(0)));
    CHECK(std::isfinite(me(1)));
}

TEST_CASE("particle filter extraction_method config switches state() output")
{
    pf_linear_dynamics dyn;
    pf_position_measurement meas;

    // Default: weighted_mean
    {
        auto pf = make_particle_filter<50>(dyn, meas, pf_config<double, 2, 1, 1>{.extraction = extraction_method::weighted_mean}, std::mt19937_64{11});

        Vector<double, 1> u = Vector<double, 1>::Zero();
        pf.predict(u);
        Vector<double, 1> z;
        z << 1.0;
        pf.update(z);

        auto st = pf.state();
        auto wm = pf.weighted_mean();
        CHECK_THAT(st(0), Catch::Matchers::WithinAbs(wm(0), 1e-14));
        CHECK_THAT(st(1), Catch::Matchers::WithinAbs(wm(1), 1e-14));
    }

    // Map mode
    {
        auto pf = make_particle_filter<50>(dyn, meas, pf_config<double, 2, 1, 1>{.extraction = extraction_method::map}, std::mt19937_64{11});

        Vector<double, 1> u = Vector<double, 1>::Zero();
        pf.predict(u);
        Vector<double, 1> z;
        z << 1.0;
        pf.update(z);

        auto st = pf.state();
        auto me = pf.map_estimate();
        CHECK_THAT(st(0), Catch::Matchers::WithinAbs(me(0), 1e-14));
        CHECK_THAT(st(1), Catch::Matchers::WithinAbs(me(1), 1e-14));
    }
}

TEST_CASE("particle filter shares dynamics_model with EKF")
{
    // The same dynamics callable satisfies dynamics_model for both EKF and PF
    auto shared_dynamics = [](const Vector<double, 2>& x, const Vector<double, 1>& u) -> Vector<double, 2>
    {
        constexpr double dt = 0.1;
        Vector<double, 2> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    };

    auto meas_fn = [](const Vector<double, 2>& x) -> Vector<double, 1>
    {
        Vector<double, 1> z;
        z(0) = x(0);
        return z;
    };

    static_assert(dynamics_model<decltype(shared_dynamics), double, 2, 1>);

    particle_filter<double, 2, 1, 1, 50, decltype(shared_dynamics), decltype(meas_fn)> pf(shared_dynamics, meas_fn, pf_config<double, 2, 1, 1>{}, systematic_resampling{}, std::mt19937_64{42});

    Vector<double, 1> u = Vector<double, 1>::Zero();
    pf.predict(u);

    Vector<double, 1> z;
    z << 1.0;
    pf.update(z);

    CHECK(std::isfinite(pf.state()(0)));
}
