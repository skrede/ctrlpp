// fixed-size templated inputs, using the belt-and-suspenders harness from
// nomalloc_harness.h: a throwing eigen_assert that survives -DNDEBUG plus a
// global allocation counter that catches heap traffic outside Eigen's own
// bookkeeping. The harness header must stay the first include of this file.
//
// Coverage: kalman_filter, ekf, ukf, mekf, manifold_ukf, complementary_filter,
// and particle_filter predict/update. Construction (including the create
// filters) happens outside the armed window; only the steady-state predict and
// update loop is guarded. The particle_filter case additionally forces the
// resampling path every step, seeds its RNG deterministically, and asserts that
// two identically seeded filters fed identical measurements stay bitwise equal.

#include "nomalloc_harness.h"
#include "hardening_helpers.h"

#include "ctrlpp/model/state_space.h"

#include "ctrlpp/estimation/ekf.h"
#include "ctrlpp/estimation/ukf.h"
#include "ctrlpp/estimation/mekf.h"
#include "ctrlpp/estimation/kalman.h"
#include "ctrlpp/estimation/manifold_ukf.h"
#include "ctrlpp/estimation/particle_filter.h"
#include "ctrlpp/estimation/complementary_filter.h"
#include "ctrlpp/estimation/estimation_types.h"
#include "ctrlpp/estimation/resampling/systematic_resampling.h"

#include "ctrlpp/lie/so3.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>

#include <random>
#include <cstddef>
#include <utility>


namespace
{

using ctrlpp::Vector;
using ctrlpp::Matrix;

template <typename Fn>
std::size_t guarded_allocations(Fn&& fn)
{
    ctrlpp_test::scoped_no_malloc guard;
    std::forward<Fn>(fn)();
    return guard.allocations();
}

struct linear_dynamics
{
    double dt = 0.1;

    auto operator()(const Vector<double, 2>& x, const Vector<double, 1>& u) const -> Vector<double, 2>
    {
        Vector<double, 2> x_next;
        x_next(0) = x(0) + dt * x(1) + 0.5 * dt * dt * u(0);
        x_next(1) = x(1) + dt * u(0);
        return x_next;
    }

    auto jacobian_x(const Vector<double, 2>&, const Vector<double, 1>&) const -> Matrix<double, 2, 2>
    {
        Matrix<double, 2, 2> f;
        f << 1.0, dt, 0.0, 1.0;
        return f;
    }

    auto jacobian_u(const Vector<double, 2>&, const Vector<double, 1>&) const -> Matrix<double, 2, 1>
    {
        Matrix<double, 2, 1> g;
        g << 0.5 * dt * dt, dt;
        return g;
    }
};

struct position_measurement
{
    auto operator()(const Vector<double, 2>& x) const -> Vector<double, 1>
    {
        return (Vector<double, 1>() << x(0)).finished();
    }

    auto jacobian(const Vector<double, 2>&) const -> Matrix<double, 1, 2>
    {
        return (Matrix<double, 1, 2>() << 1.0, 0.0).finished();
    }
};

struct gravity_measurement
{
    auto operator()(const Eigen::Quaternion<double>& q, const Vector<double, 3>&) const -> Vector<double, 3>
    {
        return q.toRotationMatrix().transpose().col(2);
    }
};

struct rotation_dynamics
{
    double dt = 0.01;

    auto operator()(const Eigen::Quaternion<double>& q, const Vector<double, 3>& omega) const
        -> Eigen::Quaternion<double>
    {
        const Vector<double, 3> phi = (omega * dt).eval();
        return (q * ctrlpp::so3::exp(phi)).normalized();
    }
};

struct attitude_measurement
{
    auto operator()(const Eigen::Quaternion<double>& q) const -> Vector<double, 3>
    {
        return q.toRotationMatrix().transpose().col(2);
    }
};

}


TEST_CASE("kalman_filter predict/update performs zero heap allocation",
          "[kalman][hardening][nomalloc]")
{
    ctrlpp::discrete_state_space<double, 2, 1, 1> sys;
    sys.A << 1.0, 0.1, 0.0, 1.0;
    sys.B << 0.005, 0.1;
    sys.C << 1.0, 0.0;
    sys.D << 0.0;

    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.01;
    Matrix<double, 1, 1> R;
    R << 1.0;

    auto kf = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(
        sys, {.Q = Q, .R = R, .x0 = Vector<double, 2>::Zero(), .P0 = Matrix<double, 2, 2>::Identity() * 10.0}));

    const Vector<double, 1> u = Vector<double, 1>::Zero();
    const Vector<double, 1> z = (Vector<double, 1>() << 0.1).finished();

    kf.predict(u);
    REQUIRE(kf.update(z).has_value());

    bool all_stepped = true;
    std::size_t allocations = 0;
    allocations = guarded_allocations([&] {
        for(int i = 0; i < 128; ++i)
        {
            kf.predict(u);
            all_stepped = all_stepped && kf.update(z).has_value();
        }
    });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(all_stepped);

    REQUIRE(allocations == 0);
    REQUIRE(kf.state().allFinite());
}

TEST_CASE("ekf predict/update performs zero heap allocation",
          "[ekf][hardening][nomalloc]")
{
    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.01;
    Matrix<double, 1, 1> R;
    R << 1.0;

    auto filter = ctrlpp::test::constructed(ctrlpp::ekf<double, 2, 1, 1, linear_dynamics, position_measurement>::create(
        linear_dynamics{}, position_measurement{},
        ctrlpp::ekf_config<double, 2, 1, 1>{
            .Q = Q, .R = R, .x0 = Vector<double, 2>::Zero(), .P0 = Matrix<double, 2, 2>::Identity() * 10.0}));

    const Vector<double, 1> u = Vector<double, 1>::Zero();
    const Vector<double, 1> z = (Vector<double, 1>() << 0.1).finished();

    filter.predict(u);
    REQUIRE(filter.update(z).has_value());

    bool all_stepped = true;
    std::size_t allocations = 0;
    allocations = guarded_allocations([&] {
        for(int i = 0; i < 128; ++i)
        {
            filter.predict(u);
            all_stepped = all_stepped && filter.update(z).has_value();
        }
    });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(all_stepped);

    REQUIRE(allocations == 0);
    REQUIRE(filter.state().allFinite());
}

TEST_CASE("ukf predict/update performs zero heap allocation",
          "[ukf][hardening][nomalloc]")
{
    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.01;
    Matrix<double, 1, 1> R;
    R << 1.0;

    auto filter = ctrlpp::test::constructed(ctrlpp::ukf<double, 2, 1, 1, linear_dynamics, position_measurement>::create(
        linear_dynamics{}, position_measurement{},
        ctrlpp::ukf_config<double, 2, 1, 1>{
            .Q = Q, .R = R, .x0 = Vector<double, 2>::Zero(), .P0 = Matrix<double, 2, 2>::Identity() * 10.0}));

    const Vector<double, 1> u = Vector<double, 1>::Zero();
    const Vector<double, 1> z = (Vector<double, 1>() << 0.1).finished();

    filter.predict(u);
    REQUIRE(filter.update(z).has_value());

    bool all_stepped = true;
    std::size_t allocations = 0;
    allocations = guarded_allocations([&] {
        for(int i = 0; i < 128; ++i)
        {
            filter.predict(u);
            all_stepped = all_stepped && filter.update(z).has_value();
        }
    });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(all_stepped);

    REQUIRE(allocations == 0);
    REQUIRE(filter.state().allFinite());
}

TEST_CASE("mekf predict/update performs zero heap allocation",
          "[mekf][hardening][nomalloc]")
{
    using mekf_type = ctrlpp::mekf<double, 3, 3, gravity_measurement>;

    ctrlpp::mekf_config<double, 3, 3> cfg;
    cfg.Q *= 1e-6;

    auto created = mekf_type::create(gravity_measurement{}, cfg);
    REQUIRE(created.has_value());
    auto& filter = *created;

    const Vector<double, 3> omega = (Vector<double, 3>() << 0.01, -0.02, 0.03).finished();
    const Vector<double, 3> z = (Vector<double, 3>() << 0.0, 0.0, 1.0).finished();

    filter.predict(omega);
    REQUIRE(filter.update(z).has_value());

    bool all_stepped = true;
    std::size_t allocations = 0;
    allocations = guarded_allocations([&] {
        for(int i = 0; i < 128; ++i)
        {
            filter.predict(omega);
            all_stepped = all_stepped && filter.update(z).has_value();
        }
    });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(all_stepped);

    REQUIRE(allocations == 0);
    REQUIRE(filter.state().allFinite());
}

TEST_CASE("manifold_ukf predict/update performs zero heap allocation",
          "[manifold_ukf][hardening][nomalloc]")
{
    using mukf_type = ctrlpp::manifold_ukf<double, 3, rotation_dynamics, attitude_measurement>;

    ctrlpp::manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-6;

    auto created = mukf_type::create(rotation_dynamics{}, attitude_measurement{}, cfg);
    REQUIRE(created.has_value());
    auto& filter = *created;

    const Vector<double, 3> omega = (Vector<double, 3>() << 0.01, -0.02, 0.03).finished();
    const Vector<double, 3> z = (Vector<double, 3>() << 0.0, 0.0, 1.0).finished();

    filter.predict(omega);
    REQUIRE(filter.update(z).has_value());

    bool all_stepped = true;
    std::size_t allocations = 0;
    allocations = guarded_allocations([&] {
        for(int i = 0; i < 128; ++i)
        {
            filter.predict(omega);
            all_stepped = all_stepped && filter.update(z).has_value();
        }
    });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(all_stepped);

    REQUIRE(allocations == 0);
    REQUIRE(filter.state().allFinite());
}

TEST_CASE("complementary_filter update performs zero heap allocation",
          "[complementary_filter][hardening][nomalloc]")
{
    ctrlpp::cf_config<double> cfg{.k_p = 2.0, .k_i = 0.005, .dt = 0.01};

    auto created = ctrlpp::complementary_filter<double>::create(cfg);
    REQUIRE(created.has_value());
    auto& filter = *created;

    const Vector<double, 3> gyro = (Vector<double, 3>() << 0.01, -0.02, 0.03).finished();
    const Vector<double, 3> accel = (Vector<double, 3>() << 0.0, 0.0, 9.81).finished();

    REQUIRE(filter.update(gyro, accel, 0.01).has_value());

    bool all_stepped = true;
    std::size_t allocations = 0;
    allocations = guarded_allocations([&] {
        for(int i = 0; i < 128; ++i)
            all_stepped = all_stepped && filter.update(gyro, accel, 0.01).has_value();
    });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(all_stepped);

    REQUIRE(allocations == 0);
    REQUIRE(filter.state().allFinite());
}

TEST_CASE("particle_filter predict/update performs zero heap allocation and is seed-deterministic",
          "[particle_filter][hardening][nomalloc]")
{
    constexpr std::size_t np = 256;

    Matrix<double, 2, 2> Q = Matrix<double, 2, 2>::Identity() * 0.01;
    Matrix<double, 1, 1> R;
    R << 0.5;

    // Force the resampling path on every update so the guard covers resample +
    // roughening (a positive ess_threshold up to NP always trips once weights
    // stop being uniform), not just the propagate/reweight step.
    const auto make_config = [&] {
        ctrlpp::pf_config<double, 2, 1, 1> cfg{};
        cfg.Q = Q;
        cfg.R = R;
        cfg.x0 = Vector<double, 2>::Zero();
        cfg.P0 = Matrix<double, 2, 2>::Identity();
        cfg.ess_threshold = static_cast<double>(np);
        return cfg;
    };

    auto pf = ctrlpp::make_particle_filter<np>(
        linear_dynamics{}, position_measurement{}, make_config(), std::mt19937_64{42});

    // Second identically seeded filter for the repeatability assertion.
    auto pf_repeat = ctrlpp::make_particle_filter<np>(
        linear_dynamics{}, position_measurement{}, make_config(), std::mt19937_64{42});

    const Vector<double, 1> u = Vector<double, 1>::Zero();
    const Vector<double, 1> z = (Vector<double, 1>() << 0.1).finished();

    pf.predict(u);
    pf.update(z);
    pf_repeat.predict(u);
    pf_repeat.update(z);

    std::size_t allocations = 0;
    allocations = guarded_allocations([&] {
        for(int i = 0; i < 64; ++i)
        {
            pf.predict(u);
            pf.update(z);
        }
    });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());

    REQUIRE(allocations == 0);

    // Repeatability: drive the twin through the identical schedule and require
    // bitwise-identical state estimates at every step.
    for(int i = 0; i < 64; ++i)
    {
        pf_repeat.predict(u);
        pf_repeat.update(z);
    }

    const Vector<double, 2> a = pf.state();
    const Vector<double, 2> b = pf_repeat.state();
    REQUIRE(a(0) == b(0));
    REQUIRE(a(1) == b(1));
}
