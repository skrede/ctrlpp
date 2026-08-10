// Compile-and-smoke tier: constructs PID, LQR, KF, EKF, UKF, MEKF, the particle
// filter, the manifold UKF, the complementary filter, SO3 and the double-S
// trajectory at Scalar=float, and drives one trivial step or evaluation on
// each. It deliberately avoids the two float-fatal absolute tolerances
// (place.h's conjugate-pair check and the biquad DC-gain guard), which the
// float runtime tier exercises separately.
//
// Every case CONSTRUCTS its type rather than naming it. An explicit class
// template instantiation instantiates member function definitions but not the
// default member initializers of a configuration aggregate nothing builds, and
// those initializers are where the single-precision conversion diagnostics sit,
// so a tier that only named its types would read clean over code no compiler
// had looked at.

#include "hardening_helpers.h"
#include "ctrlpp/lie/so3.h"
#include "ctrlpp/control/lqr.h"
#include "ctrlpp/control/pid.h"
#include "ctrlpp/estimation/ekf.h"
#include "ctrlpp/estimation/ukf.h"
#include "ctrlpp/estimation/mekf.h"
#include "ctrlpp/estimation/kalman.h"
#include "ctrlpp/estimation/manifold_ukf.h"
#include "ctrlpp/estimation/particle_filter.h"
#include "ctrlpp/estimation/complementary_filter.h"
#include "ctrlpp/trajectory/double_s_trajectory.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <random>
#include <cstddef>

using namespace ctrlpp;

namespace
{

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr std::size_t NY = 1;

struct linear_dynamics_f
{
    auto operator()(const Vector<float, NX>& x, const Vector<float, NU>& u) const -> Vector<float, NX>
    {
        Vector<float, NX> x_next;
        x_next(0) = x(0) + 0.1f * x(1);
        x_next(1) = x(1) + 0.1f * u(0);
        return x_next;
    }
};

struct position_measurement_f
{
    auto operator()(const Vector<float, NX>& x) const -> Vector<float, NY>
    {
        return (Vector<float, NY>() << x(0)).finished();
    }
};

}

TEST_CASE("PID instantiates and steps at Scalar=float", "[float][anchor]")
{
    pid_config<float, 1> cfg;
    cfg.kp = Vector<float, 1>::Constant(1.0f);
    cfg.ki = Vector<float, 1>::Constant(0.1f);
    cfg.kd = Vector<float, 1>::Constant(0.01f);

    pid<float, 1> controller(cfg);
    auto u = ctrlpp::test::commanded(controller.compute(Vector<float, 1>::Constant(1.0f), Vector<float, 1>::Zero(), 0.1f));

    REQUIRE(std::isfinite(u(0)));
}

TEST_CASE("LQR instantiates and computes at Scalar=float", "[float][anchor]")
{
    // Constructed directly from a precomputed gain (no DARE solve involved):
    // the LQR class itself carries no float-fatal tolerance.
    Eigen::Matrix<float, 1, 2> K;
    K << 1.0f, 0.5f;
    lqr<float, 2, 1> controller(K);

    Vector<float, 2> x;
    x << 1.0f, 0.0f;
    auto u = controller.compute(x);

    REQUIRE(std::isfinite(u(0)));
}

TEST_CASE("KF instantiates and steps at Scalar=float", "[float][anchor]")
{
    discrete_state_space<float, NX, NU, NY> sys;
    sys.A << 1.0f, 0.1f, 0.0f, 1.0f;
    sys.B << 0.005f, 0.1f;
    sys.C << 1.0f, 0.0f;
    sys.D = Matrix<float, NY, NU>::Zero();

    kalman_config<float, NX, NU, NY> cfg;
    auto filt = ctrlpp::test::constructed(kalman_filter<float, NX, NU, NY>::create(sys, cfg));

    filt.predict(Vector<float, NU>::Zero());
    REQUIRE(filt.update(Vector<float, NY>::Constant(0.5f)).has_value());

    REQUIRE(filt.state().allFinite());
    REQUIRE(filt.covariance().allFinite());
}

TEST_CASE("EKF instantiates and steps at Scalar=float", "[float][anchor]")
{
    ekf_config<float, NX, NU, NY> cfg;
    auto filt = ctrlpp::test::constructed(ekf<float, NX, NU, NY, linear_dynamics_f, position_measurement_f>::create(linear_dynamics_f{}, position_measurement_f{}, cfg));

    filt.predict(Vector<float, NU>::Zero());
    REQUIRE(filt.update(Vector<float, NY>::Constant(0.5f)).has_value());

    REQUIRE(filt.state().allFinite());
    REQUIRE(filt.covariance().allFinite());
}

TEST_CASE("UKF instantiates and steps at Scalar=float", "[float][anchor]")
{
    ukf_config<float, NX, NU, NY> cfg;
    auto filt = ctrlpp::test::constructed(ukf<float, NX, NU, NY, linear_dynamics_f, position_measurement_f>::create(linear_dynamics_f{}, position_measurement_f{}, cfg));

    filt.predict(Vector<float, NU>::Zero());
    REQUIRE(filt.update(Vector<float, NY>::Constant(0.5f)).has_value());

    REQUIRE(filt.state().allFinite());
    REQUIRE(filt.covariance().allFinite());
}

TEST_CASE("MEKF instantiates and steps at Scalar=float", "[float][anchor]")
{
    constexpr std::size_t NB = 3;
    constexpr std::size_t NY_MEKF = 3;

    struct vector_observation_f
    {
        auto operator()(const Eigen::Quaternion<float>& q, const Vector<float, NB>& /*b*/) const -> Vector<float, NY_MEKF>
        {
            Vector<float, 3> r_w{0.0f, 0.0f, 1.0f};
            return q.conjugate() * r_w;
        }
    };

    mekf_config<float, NB, NY_MEKF> cfg;
    cfg.dt = 0.1f;
    auto filt_result = mekf<float, NB, NY_MEKF, vector_observation_f>::create(vector_observation_f{}, cfg);
    REQUIRE(filt_result.has_value());
    auto& filt = *filt_result;

    Vector<float, 3> omega{0.1f, 0.0f, 0.0f};
    filt.predict(omega, 0.1f);
    REQUIRE(filt.update(Vector<float, NY_MEKF>{0.0f, 0.0f, 1.0f}).has_value());

    REQUIRE(filt.state().allFinite());
    REQUIRE(filt.covariance().allFinite());
}

TEST_CASE("particle filter instantiates and steps at Scalar=float", "[float][anchor]")
{
    constexpr std::size_t NP = 16;

    auto filt = make_particle_filter<NP>(
        linear_dynamics_f{}, position_measurement_f{}, pf_config<float, NX, NU, NY>{}, std::mt19937_64{42U});

    filt.predict(Vector<float, NU>::Zero());
    filt.update(Vector<float, NY>::Constant(0.5f));

    REQUIRE(filt.state().allFinite());
}

TEST_CASE("manifold UKF instantiates and steps at Scalar=float", "[float][anchor]")
{
    constexpr std::size_t NY_MUKF = 3;

    struct attitude_propagation_f
    {
        auto operator()(const Eigen::Quaternion<float>& q, const Vector<float, 3>& /*omega*/) const -> Eigen::Quaternion<float>
        {
            return q;
        }
    };

    struct gravity_observation_f
    {
        auto operator()(const Eigen::Quaternion<float>& q) const -> Vector<float, NY_MUKF>
        {
            Vector<float, 3> r_w{0.0f, 0.0f, 1.0f};
            return q.conjugate() * r_w;
        }
    };

    auto filt_result = manifold_ukf<float, NY_MUKF, attitude_propagation_f, gravity_observation_f>::create(
        attitude_propagation_f{}, gravity_observation_f{}, manifold_ukf_config<float, NY_MUKF>{});
    REQUIRE(filt_result.has_value());
    auto& filt = *filt_result;

    filt.predict(Vector<float, 3>{0.1f, 0.0f, 0.0f});
    REQUIRE(filt.update(Vector<float, NY_MUKF>{0.0f, 0.0f, 1.0f}).has_value());

    REQUIRE(filt.state().allFinite());
    REQUIRE(filt.covariance().allFinite());
}

TEST_CASE("complementary_filter instantiates and steps at Scalar=float", "[float][anchor]")
{
    cf_config<float> cfg;
    auto filt_result = complementary_filter<float>::create(cfg);
    REQUIRE(filt_result.has_value());
    auto& filt = *filt_result;

    Vector<float, 3> gyro{0.1f, 0.0f, 0.0f};
    Vector<float, 3> accel{0.0f, 0.0f, 9.8f};
    REQUIRE(filt.update(gyro, accel, 0.01f).has_value());

    REQUIRE(filt.state().allFinite());
}

TEST_CASE("SO3 exp/log round-trip at Scalar=float", "[float][anchor]")
{
    Vector<float, 3> phi{0.2f, -0.1f, 0.05f};
    auto q = so3::exp(phi);
    auto phi_back = so3::log(q);

    REQUIRE(phi_back.allFinite());
    REQUIRE(q.coeffs().allFinite());
}

TEST_CASE("double-S trajectory instantiates and evaluates at Scalar=float", "[float][anchor]")
{
    auto traj_result = double_s_trajectory<float>::create(
        {.q0 = 0.0f, .q1 = 1.0f, .v_max = 1.0f, .a_max = 1.0f, .j_max = 1.0f});
    REQUIRE(traj_result.has_value());
    auto& traj = *traj_result;

    const auto point = traj.evaluate(traj.duration() * 0.5f);

    REQUIRE(point.position.allFinite());
    REQUIRE(point.velocity.allFinite());
    REQUIRE(point.acceleration.allFinite());
}
