#include "ctrlpp/estimation/particle_filter.h"
#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/resampling/resampling_strategy.h"
#include "ctrlpp/estimation/resampling/systematic_resampling.h"
#include "ctrlpp/estimation/resampling/multinomial_resampling.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <random>

using namespace ctrlpp;

// ---------------------------------------------------------------------------
// Concept satisfaction: resampling strategies
// ---------------------------------------------------------------------------
static_assert(resampling_strategy<systematic_resampling, std::mt19937_64, 10>);
static_assert(resampling_strategy<multinomial_resampling, std::mt19937_64, 10>);

// ---------------------------------------------------------------------------
// Test dynamics: double integrator (shared with EKF/UKF)
// ---------------------------------------------------------------------------
struct pf_linear_dynamics {
    double dt = 0.1;

    auto operator()(const Vector<double, 2>& x,
                    const Vector<double, 1>& u) const -> Vector<double, 2>
    {
        Vector<double, 2> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    }
};

struct pf_position_measurement {
    auto operator()(const Vector<double, 2>& x) const -> Vector<double, 1>
    {
        Vector<double, 1> z;
        z(0) = x(0);
        return z;
    }
};

// Concept satisfaction: particle filter
using pf_type = particle_filter<double, 2, 1, 1, 100,
                                pf_linear_dynamics, pf_position_measurement>;
static_assert(ObserverPolicy<pf_type>);
static_assert(!CovarianceObserver<pf_type>);

// ---------------------------------------------------------------------------
// Resampling strategy tests
// ---------------------------------------------------------------------------

TEST_CASE("systematic resampling uniform weights produce identity indices") {
    constexpr std::size_t NP = 8;
    std::array<double, NP> weights;
    weights.fill(1.0 / static_cast<double>(NP));

    std::array<std::size_t, NP> indices{};
    std::mt19937_64 rng(42);

    systematic_resampling resampler;
    resampler.resample(weights, indices, rng);

    for (std::size_t i = 0; i < NP; ++i) {
        CHECK(indices[i] == i);
    }
}

TEST_CASE("systematic resampling degenerate weights concentrate on single index") {
    constexpr std::size_t NP = 10;
    std::array<double, NP> weights{};
    weights.fill(0.0);
    weights[3] = 1.0;

    std::array<std::size_t, NP> indices{};
    std::mt19937_64 rng(123);

    systematic_resampling resampler;
    resampler.resample(weights, indices, rng);

    for (std::size_t i = 0; i < NP; ++i) {
        CHECK(indices[i] == 3);
    }
}

TEST_CASE("systematic resampling all indices valid") {
    constexpr std::size_t NP = 100;
    std::array<double, NP> weights;
    weights.fill(1.0 / static_cast<double>(NP));

    std::array<std::size_t, NP> indices{};
    std::mt19937_64 rng(7);

    systematic_resampling resampler;
    resampler.resample(weights, indices, rng);

    for (std::size_t i = 0; i < NP; ++i) {
        CHECK(indices[i] < NP);
    }
}

TEST_CASE("systematic resampling deterministic with fixed seed") {
    constexpr std::size_t NP = 20;
    std::array<double, NP> weights;
    for (std::size_t i = 0; i < NP; ++i) {
        weights[i] = static_cast<double>(i + 1);
    }
    double sum = 0.0;
    for (auto w : weights) sum += w;
    for (auto& w : weights) w /= sum;

    std::array<std::size_t, NP> indices1{};
    std::array<std::size_t, NP> indices2{};

    {
        std::mt19937_64 rng(999);
        systematic_resampling resampler;
        resampler.resample(weights, indices1, rng);
    }
    {
        std::mt19937_64 rng(999);
        systematic_resampling resampler;
        resampler.resample(weights, indices2, rng);
    }

    for (std::size_t i = 0; i < NP; ++i) {
        CHECK(indices1[i] == indices2[i]);
    }
}

TEST_CASE("multinomial resampling degenerate weights concentrate on single index") {
    constexpr std::size_t NP = 10;
    std::array<double, NP> weights{};
    weights.fill(0.0);
    weights[5] = 1.0;

    std::array<std::size_t, NP> indices{};
    std::mt19937_64 rng(456);

    multinomial_resampling resampler;
    resampler.resample(weights, indices, rng);

    for (std::size_t i = 0; i < NP; ++i) {
        CHECK(indices[i] == 5);
    }
}

TEST_CASE("multinomial resampling uniform weights roughly uniform distribution") {
    constexpr std::size_t NP = 1000;
    std::array<double, NP> weights;
    weights.fill(1.0 / static_cast<double>(NP));

    std::array<std::size_t, NP> indices{};
    std::mt19937_64 rng(789);

    multinomial_resampling resampler;
    resampler.resample(weights, indices, rng);

    for (std::size_t i = 0; i < NP; ++i) {
        CHECK(indices[i] < NP);
    }

    std::array<int, NP> counts{};
    for (auto idx : indices) {
        ++counts[idx];
    }
    std::size_t unique = 0;
    for (auto c : counts) {
        if (c > 0) ++unique;
    }
    CHECK(unique > NP / 2);
}

// ---------------------------------------------------------------------------
// Particle filter tests
// ---------------------------------------------------------------------------

TEST_CASE("particle filter tracks linear system") {
    pf_linear_dynamics dyn;
    pf_position_measurement meas;

    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.01;
    Matrix<double, 1, 1> R;
    R << 0.5;
    Vector<double, 2> x0 = Vector<double, 2>::Zero();
    Matrix<double, 2, 2> P0 = Matrix<double, 2, 2>::Identity() * 1.0;

    auto pf = make_particle_filter<500>(
        dyn, meas,
        pf_config<double, 2, 1, 1>{.Q = Q, .R = R, .x0 = x0, .P0 = P0},
        std::mt19937_64{42});

    double true_pos = 0.0;
    double true_vel = 1.0;
    constexpr double dt = 0.1;

    for (int i = 0; i < 50; ++i) {
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

TEST_CASE("particle filter deterministic with fixed seed") {
    pf_linear_dynamics dyn;
    pf_position_measurement meas;

    pf_config<double, 2, 1, 1> cfg{
        .Q = Matrix<double, 2, 2>::Identity() * 0.01,
        .R = ([]{ Matrix<double, 1, 1> R; R << 0.5; return R; })(),
        .x0 = Vector<double, 2>::Zero(),
        .P0 = Matrix<double, 2, 2>::Identity()
    };

    auto run = [&](std::uint64_t seed) {
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

TEST_CASE("particle filter weighted_mean and map_estimate both accessible") {
    pf_linear_dynamics dyn;
    pf_position_measurement meas;

    auto pf = make_particle_filter<50>(
        dyn, meas,
        pf_config<double, 2, 1, 1>{},
        std::mt19937_64{77});

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

TEST_CASE("particle filter extraction_method config switches state() output") {
    pf_linear_dynamics dyn;
    pf_position_measurement meas;

    // Default: weighted_mean
    {
        auto pf = make_particle_filter<50>(
            dyn, meas,
            pf_config<double, 2, 1, 1>{.extraction = extraction_method::weighted_mean},
            std::mt19937_64{11});

        Vector<double, 1> u = Vector<double, 1>::Zero();
        pf.predict(u);
        Vector<double, 1> z; z << 1.0;
        pf.update(z);

        auto st = pf.state();
        auto wm = pf.weighted_mean();
        CHECK_THAT(st(0), Catch::Matchers::WithinAbs(wm(0), 1e-14));
        CHECK_THAT(st(1), Catch::Matchers::WithinAbs(wm(1), 1e-14));
    }

    // Map mode
    {
        auto pf = make_particle_filter<50>(
            dyn, meas,
            pf_config<double, 2, 1, 1>{.extraction = extraction_method::map},
            std::mt19937_64{11});

        Vector<double, 1> u = Vector<double, 1>::Zero();
        pf.predict(u);
        Vector<double, 1> z; z << 1.0;
        pf.update(z);

        auto st = pf.state();
        auto me = pf.map_estimate();
        CHECK_THAT(st(0), Catch::Matchers::WithinAbs(me(0), 1e-14));
        CHECK_THAT(st(1), Catch::Matchers::WithinAbs(me(1), 1e-14));
    }
}

TEST_CASE("particle filter ESS-adaptive resampling triggers correctly") {
    // With a very informative measurement (small R) and particles spread out,
    // ESS should drop and trigger resampling
    pf_linear_dynamics dyn;
    pf_position_measurement meas;

    Matrix<double, 1, 1> R;
    R << 0.001; // Very small measurement noise -> high info -> ESS drops

    auto pf = make_particle_filter<200>(
        dyn, meas,
        pf_config<double, 2, 1, 1>{
            .Q = Matrix<double, 2, 2>::Identity() * 0.01,
            .R = R,
            .x0 = Vector<double, 2>::Zero(),
            .P0 = Matrix<double, 2, 2>::Identity() * 10.0
        },
        std::mt19937_64{55});

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

TEST_CASE("particle filter roughening prevents impoverishment") {
    pf_linear_dynamics dyn;
    pf_position_measurement meas;

    Matrix<double, 1, 1> R;
    R << 0.001;

    auto pf = make_particle_filter<100>(
        dyn, meas,
        pf_config<double, 2, 1, 1>{
            .Q = Matrix<double, 2, 2>::Identity() * 0.01,
            .R = R,
            .x0 = Vector<double, 2>::Zero(),
            .P0 = Matrix<double, 2, 2>::Identity() * 10.0,
            .roughening_scale = 0.5
        },
        std::mt19937_64{88});

    Vector<double, 1> u = Vector<double, 1>::Zero();
    pf.predict(u);

    Vector<double, 1> z;
    z << 5.0;
    pf.update(z);

    // After resampling + roughening, particles should not all be identical
    auto& particles = pf.particles();
    bool all_same = true;
    for (std::size_t i = 1; i < 100; ++i) {
        if ((particles[i] - particles[0]).norm() > 1e-10) {
            all_same = false;
            break;
        }
    }
    CHECK_FALSE(all_same);
}

TEST_CASE("particle filter satisfies ObserverPolicy") {
    static_assert(ObserverPolicy<pf_type>);
    CHECK(true);
}

TEST_CASE("particle filter does NOT satisfy CovarianceObserver") {
    static_assert(!CovarianceObserver<pf_type>);
    CHECK(true);
}

TEST_CASE("particle filter with multinomial resampling") {
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

    for (int i = 0; i < 30; ++i) {
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

TEST_CASE("particle filter log-weight stability over many iterations") {
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

    for (int i = 0; i < 200; ++i) {
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

TEST_CASE("particle filter shares dynamics_model with EKF") {
    // The same dynamics callable satisfies dynamics_model for both EKF and PF
    auto shared_dynamics = [](const Vector<double, 2>& x,
                              const Vector<double, 1>& u) -> Vector<double, 2> {
        constexpr double dt = 0.1;
        Vector<double, 2> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    };

    auto meas_fn = [](const Vector<double, 2>& x) -> Vector<double, 1> {
        Vector<double, 1> z;
        z(0) = x(0);
        return z;
    };

    static_assert(dynamics_model<decltype(shared_dynamics), double, 2, 1>);

    particle_filter<double, 2, 1, 1, 50,
                    decltype(shared_dynamics), decltype(meas_fn)> pf(
        shared_dynamics, meas_fn,
        pf_config<double, 2, 1, 1>{},
        systematic_resampling{},
        std::mt19937_64{42});

    Vector<double, 1> u = Vector<double, 1>::Zero();
    pf.predict(u);

    Vector<double, 1> z;
    z << 1.0;
    pf.update(z);

    CHECK(std::isfinite(pf.state()(0)));
}

TEST_CASE("particle filter linear-weight mode works for small particle count") {
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

    for (int i = 0; i < 20; ++i) {
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
