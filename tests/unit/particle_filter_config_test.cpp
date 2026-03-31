#include "ctrlpp/estimation/particle_filter.h"
#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/resampling/multinomial_resampling.h"
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

// Concept satisfaction: particle filter
using pf_type = particle_filter<double, 2, 1, 1, 100, pf_linear_dynamics, pf_position_measurement>;
static_assert(ObserverPolicy<pf_type>);
static_assert(!CovarianceObserver<pf_type>);

// ---------------------------------------------------------------------------
// Configuration and policy tests
// ---------------------------------------------------------------------------

TEST_CASE("particle filter ESS-adaptive resampling triggers correctly")
{
    // With a very informative measurement (small R) and particles spread out,
    // ESS should drop and trigger resampling
    pf_linear_dynamics dyn;
    pf_position_measurement meas;

    Matrix<double, 1, 1> R;
    R << 0.001; // Very small measurement noise -> high info -> ESS drops

    auto pf = make_particle_filter<200>(
        dyn, meas, pf_config<double, 2, 1, 1>{.Q = Matrix<double, 2, 2>::Identity() * 0.01, .R = R, .x0 = Vector<double, 2>::Zero(), .P0 = Matrix<double, 2, 2>::Identity() * 10.0}, std::mt19937_64{55});

    Vector<double, 1> u = Vector<double, 1>::Zero();
    pf.predict(u);

    Vector<double, 1> z;
    z << 5.0;
    pf.update(z);

    // After update with highly informative measurement, the filter should still
    // produce a finite estimate (resampling should have occurred)
    auto est = pf.state();
    CHECK(std::isfinite(est(0)));
    CHECK(std::isfinite(est(1)));
}

TEST_CASE("particle filter roughening prevents impoverishment")
{
    pf_linear_dynamics dyn;
    pf_position_measurement meas;

    Matrix<double, 1, 1> R;
    R << 0.001;

    auto pf =
        make_particle_filter<100>(dyn,
                                  meas,
                                  pf_config<double, 2, 1, 1>{.Q = Matrix<double, 2, 2>::Identity() * 0.01, .R = R, .x0 = Vector<double, 2>::Zero(), .P0 = Matrix<double, 2, 2>::Identity() * 10.0, .roughening_scale = 0.5},
                                  std::mt19937_64{88});

    Vector<double, 1> u = Vector<double, 1>::Zero();
    pf.predict(u);

    Vector<double, 1> z;
    z << 5.0;
    pf.update(z);

    // After resampling + roughening, particles should not all be identical
    auto& particles = pf.particles();
    bool all_same = true;
    for(std::size_t i = 1; i < 100; ++i)
    {
        if((particles[i] - particles[0]).norm() > 1e-10)
        {
            all_same = false;
            break;
        }
    }
    CHECK_FALSE(all_same);
}

TEST_CASE("particle filter satisfies ObserverPolicy")
{
    static_assert(ObserverPolicy<pf_type>);
    CHECK(true);
}

TEST_CASE("particle filter does NOT satisfy CovarianceObserver")
{
    static_assert(!CovarianceObserver<pf_type>);
    CHECK(true);
}

TEST_CASE("particle filter with multinomial resampling")
{
    pf_linear_dynamics dyn;
    pf_position_measurement meas;

    particle_filter<double, 2, 1, 1, 300,
                    pf_linear_dynamics, pf_position_measurement,
                    multinomial_resampling> pf(
        dyn, meas,
        pf_config<double, 2, 1, 1>{
            .Q = Matrix<double, 2, 2>::Identity() * 0.01,
            .R = ([]{ Matrix<double, 1, 1> R; R << 0.5; return R; })(),
            .x0 = Vector<double, 2>::Zero(),
            .P0 = Matrix<double, 2, 2>::Identity()
        },
        multinomial_resampling{},
        std::mt19937_64{321});

    double true_pos = 0.0;
    double true_vel = 1.0;
    constexpr double dt = 0.1;

    for(int i = 0; i < 30; ++i)
    {
        true_pos += true_vel * dt;

        Vector<double, 1> u = Vector<double, 1>::Zero();
        pf.predict(u);

        Vector<double, 1> z;
        z << true_pos + 0.1 * std::sin(static_cast<double>(i));
        pf.update(z);
    }

    auto est = pf.state();
    CHECK_THAT(est(0), Catch::Matchers::WithinAbs(true_pos, 1.5));
}

TEST_CASE("particle filter log-weight stability over many iterations")
{
    pf_linear_dynamics dyn;
    pf_position_measurement meas;

    auto pf = make_particle_filter<200>(
        dyn, meas,
        pf_config<double, 2, 1, 1>{
            .Q = Matrix<double, 2, 2>::Identity() * 0.01,
            .R = ([]{ Matrix<double, 1, 1> R; R << 0.5; return R; })(),
            .weights = weight_representation::log
        },
        std::mt19937_64{999});

    for(int i = 0; i < 200; ++i)
    {
        Vector<double, 1> u = Vector<double, 1>::Zero();
        pf.predict(u);

        Vector<double, 1> z;
        z << static_cast<double>(i) * 0.1;
        pf.update(z);

        auto est = pf.state();
        CHECK(std::isfinite(est(0)));
        CHECK(std::isfinite(est(1)));
    }
}

TEST_CASE("particle filter linear-weight mode works for small particle count")
{
    pf_linear_dynamics dyn;
    pf_position_measurement meas;

    auto pf = make_particle_filter<50>(
        dyn, meas,
        pf_config<double, 2, 1, 1>{
            .Q = Matrix<double, 2, 2>::Identity() * 0.01,
            .R = ([]{ Matrix<double, 1, 1> R; R << 0.5; return R; })(),
            .weights = weight_representation::linear
        },
        std::mt19937_64{777});

    for(int i = 0; i < 20; ++i)
    {
        Vector<double, 1> u = Vector<double, 1>::Zero();
        pf.predict(u);

        Vector<double, 1> z;
        z << static_cast<double>(i) * 0.1;
        pf.update(z);

        auto est = pf.state();
        CHECK(std::isfinite(est(0)));
        CHECK(std::isfinite(est(1)));
    }
}
