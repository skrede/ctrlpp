#include "ctrlpp/estimation/manifold_ukf.h"
#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/estimation/estimation_types.h"
#include "ctrlpp/lie/so3.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <utility>

using namespace ctrlpp;
using Catch::Matchers::WithinAbs;

namespace
{

struct simple_rotation_dynamics
{
    double dt = 0.01;

    auto operator()(const Eigen::Quaternion<double>& q, const Vector<double, 3>& omega) const -> Eigen::Quaternion<double>
    {
        Vector<double, 3> phi = (omega * dt).eval();
        return (q * so3::exp(phi)).normalized();
    }
};

struct gravity_meas
{
    auto operator()(const Eigen::Quaternion<double>& q) const -> Vector<double, 3>
    {
        return q.toRotationMatrix().transpose().col(2);
    }
};

using MukfType = manifold_ukf<double, 3, simple_rotation_dynamics, gravity_meas>;

// create() is the only construction path and it is fallible, so every
// valid-input site goes through it and asserts success here. The rejection
// cases below do not use this helper: they assert the specific enumerator.
auto make_filter(simple_rotation_dynamics dyn, const manifold_ukf_config<double, 3>& cfg) -> MukfType
{
    auto created = MukfType::create(std::move(dyn), gravity_meas{}, cfg);
    REQUIRE(created.has_value());
    return *std::move(created);
}

} // namespace

TEST_CASE("manifold_ukf predict with zero angular velocity preserves attitude", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-6;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    for(int i = 0; i < 100; ++i)
        filter.predict(Vector<double, 3>::Zero());

    auto q = filter.attitude();
    double angle = 2.0 * std::acos(std::clamp(std::abs(q.w()), 0.0, 1.0));
    REQUIRE(angle < 0.01);
}

TEST_CASE("manifold_ukf update corrects attitude toward measurement", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.q0 = so3::exp(Vector<double, 3>{0.174, 0.0, 0.0});
    cfg.Q *= 1e-6;
    cfg.R *= 0.01;
    cfg.P0 *= 10.0;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    Vector<double, 3> gravity_world{0.0, 0.0, 1.0};

    for(int i = 0; i < 50; ++i)
    {
        filter.predict(Vector<double, 3>::Zero());
        REQUIRE(filter.update(gravity_world).has_value());
    }

    auto q = filter.attitude();
    double angle = 2.0 * std::acos(std::clamp(std::abs(q.w()), 0.0, 1.0));
    REQUIRE(angle < 0.1);
}

TEST_CASE("manifold_ukf covariance stays symmetric and PSD", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-4;
    cfg.R *= 0.1;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    Vector<double, 3> omega{0.01, -0.02, 0.03};
    Vector<double, 3> z{0.0, 0.0, 1.0};

    for(int i = 0; i < 100; ++i)
    {
        filter.predict(omega);
        REQUIRE(filter.update(z).has_value());
    }

    auto P = filter.covariance();
    double asym = (P - P.transpose()).norm();
    REQUIRE(asym < 1e-10);

    Eigen::SelfAdjointEigenSolver<Matrix<double, 3, 3>> es(P);
    REQUIRE(es.eigenvalues().minCoeff() >= -1e-10);
}

TEST_CASE("manifold_ukf geodesic mean converges", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-6;
    cfg.R *= 0.01;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    for(int i = 0; i < 30; ++i)
    {
        filter.predict(Vector<double, 3>{0.01, 0.0, 0.0});
        REQUIRE(filter.update(Vector<double, 3>{0.0, 0.0, 1.0}).has_value());
    }

    auto q = filter.attitude();
    double qnorm = std::sqrt(q.w() * q.w() + q.x() * q.x() + q.y() * q.y() + q.z() * q.z());
    REQUIRE_THAT(qnorm, WithinAbs(1.0, 1e-10));
}

TEST_CASE("manifold_ukf innovation is finite", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    filter.predict(Vector<double, 3>{0.01, 0.0, 0.0});
    REQUIRE(filter.update(Vector<double, 3>{0.0, 0.0, 1.0}).has_value());

    REQUIRE(std::isfinite(filter.innovation().norm()));
}

TEST_CASE("manifold_ukf satisfies ObserverPolicy concept", "[manifold_ukf]")
{
    static_assert(ObserverPolicy<MukfType>);
    static_assert(CovarianceObserver<MukfType>);
}

TEST_CASE("manifold_ukf tracks constant rotation", "[manifold_ukf]")
{
    simple_rotation_dynamics dyn;
    dyn.dt = 0.01;

    manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-8;
    cfg.R *= 0.01;
    cfg.P0 *= 10.0;
    auto filter = make_filter(dyn, cfg);

    Vector<double, 3> omega{0.0, 0.0, 0.1};
    Vector<double, 3> gravity_world{0.0, 0.0, 1.0};

    // Run 100 steps: expected rotation = 100 * 0.01 * 0.1 = 0.1 rad about z
    for(int i = 0; i < 100; ++i)
    {
        filter.predict(omega);
        // Generate measurement from ground truth
        Vector<double, 3> phi_true = (omega * dyn.dt * static_cast<double>(i + 1)).eval();
        Eigen::Quaternion<double> q_true = so3::exp(phi_true);
        Vector<double, 3> z = q_true.toRotationMatrix().transpose().col(2);
        REQUIRE(filter.update(z).has_value());
    }

    auto q_est = filter.attitude();
    Vector<double, 3> phi = so3::log(q_est);
    REQUIRE(std::abs(phi(2)) > 0.01); // Some rotation accumulated
}

// ---------------------------------------------------------------------------
// Edge / failure-path tests
// ---------------------------------------------------------------------------

TEST_CASE("manifold_ukf predict-only grows covariance", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.Q = Matrix<double, 3, 3>::Identity() * 0.01;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    auto P_before = filter.covariance();
    double trace_before = P_before.trace();

    for(int i = 0; i < 50; ++i)
        filter.predict(Vector<double, 3>::Zero());

    double trace_after = filter.covariance().trace();
    // Covariance should grow when only predicting (no measurement updates)
    REQUIRE(trace_after > trace_before);
}

TEST_CASE("manifold_ukf predict-only covariance stays symmetric PSD", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.Q = Matrix<double, 3, 3>::Identity() * 0.001;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    for(int i = 0; i < 100; ++i)
        filter.predict(Vector<double, 3>{0.1, -0.05, 0.02});

    auto P = filter.covariance();
    double asym = (P - P.transpose()).norm();
    REQUIRE(asym < 1e-10);

    Eigen::SelfAdjointEigenSolver<Matrix<double, 3, 3>> es(P);
    REQUIRE(es.eigenvalues().minCoeff() >= -1e-10);
}

TEST_CASE("manifold_ukf large initial error converges with updates", "[manifold_ukf]")
{
    // Start with ~170 degree error (near antipodal)
    manifold_ukf_config<double, 3> cfg;
    cfg.q0 = so3::exp(Vector<double, 3>{2.9, 0.0, 0.0}); // ~166 degrees
    cfg.Q *= 1e-6;
    cfg.R *= 0.01;
    cfg.P0 *= 100.0; // Very uncertain
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    Vector<double, 3> gravity_world{0.0, 0.0, 1.0};

    for(int i = 0; i < 200; ++i)
    {
        filter.predict(Vector<double, 3>::Zero());
        REQUIRE(filter.update(gravity_world).has_value());
    }

    auto q = filter.attitude();
    // Should have made progress toward identity (gravity-aligned)
    double angle = 2.0 * std::acos(std::clamp(std::abs(q.w()), 0.0, 1.0));
    // With large initial error the filter should at least reduce the error significantly
    REQUIRE(angle < 2.0);
}

TEST_CASE("manifold_ukf with very small process noise", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.Q = Matrix<double, 3, 3>::Identity() * 1e-12;
    cfg.R *= 0.1;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    // The filter should still work without numerical issues
    for(int i = 0; i < 50; ++i)
    {
        filter.predict(Vector<double, 3>::Zero());
        REQUIRE(filter.update(Vector<double, 3>{0.0, 0.0, 1.0}).has_value());
    }

    auto q = filter.attitude();
    REQUIRE(std::isfinite(q.w()));
    REQUIRE(std::isfinite(q.x()));
    REQUIRE(std::isfinite(q.y()));
    REQUIRE(std::isfinite(q.z()));

    auto P = filter.covariance();
    REQUIRE(std::isfinite(P.norm()));
}

TEST_CASE("manifold_ukf with very large measurement noise trusts prediction", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-6;
    cfg.R = Matrix<double, 3, 3>::Identity() * 1e6; // Huge measurement noise
    cfg.P0 *= 0.001; // Very confident initial state
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    // Give a measurement that disagrees with identity orientation
    for(int i = 0; i < 20; ++i)
    {
        filter.predict(Vector<double, 3>::Zero());
        // Wrong direction, but a perfectly finite measurement, so the step
        // succeeds and the filter is expected to disagree with it.
        REQUIRE(filter.update(Vector<double, 3>{0.5, 0.5, 0.5}).has_value());
    }

    // Filter should mostly ignore measurement due to high R
    auto q = filter.attitude();
    double angle = 2.0 * std::acos(std::clamp(std::abs(q.w()), 0.0, 1.0));
    REQUIRE(angle < 0.5); // Still close to identity
}

TEST_CASE("manifold_ukf geodesic mean with max_iter=1", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.geodesic_mean_max_iter = 1;
    cfg.Q *= 1e-4;
    cfg.R *= 0.1;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    // Should still produce finite results even with single geodesic iteration
    for(int i = 0; i < 20; ++i)
    {
        filter.predict(Vector<double, 3>{0.01, 0.0, 0.0});
        REQUIRE(filter.update(Vector<double, 3>{0.0, 0.0, 1.0}).has_value());
    }

    auto q = filter.attitude();
    REQUIRE(std::isfinite(q.w()));
    double qnorm = std::sqrt(q.w() * q.w() + q.x() * q.x() + q.y() * q.y() + q.z() * q.z());
    REQUIRE_THAT(qnorm, WithinAbs(1.0, 1e-10));
}

TEST_CASE("manifold_ukf with loose geodesic tolerance", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.geodesic_mean_tol = 1.0; // Very loose -- should converge immediately
    cfg.Q *= 1e-4;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    for(int i = 0; i < 10; ++i)
        filter.predict(Vector<double, 3>{0.01, -0.01, 0.02});

    auto q = filter.attitude();
    REQUIRE(std::isfinite(q.w()));
}

TEST_CASE("manifold_ukf state cache matches attitude quaternion", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    filter.predict(Vector<double, 3>{0.05, -0.03, 0.01});
    REQUIRE(filter.update(Vector<double, 3>{0.0, 0.0, 1.0}).has_value());

    auto q = filter.attitude();
    auto s = filter.state();

    // state_cache should be so3::to_vec(q) = [w, x, y, z]
    REQUIRE_THAT(s(0), WithinAbs(q.w(), 1e-15));
    REQUIRE_THAT(s(1), WithinAbs(q.x(), 1e-15));
    REQUIRE_THAT(s(2), WithinAbs(q.y(), 1e-15));
    REQUIRE_THAT(s(3), WithinAbs(q.z(), 1e-15));
}

TEST_CASE("manifold_ukf innovation decreases as filter converges", "[manifold_ukf]")
{
    manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-8;
    cfg.R *= 0.01;
    auto filter = make_filter(simple_rotation_dynamics{}, cfg);

    Vector<double, 3> gravity{0.0, 0.0, 1.0};

    // Run several predict/update cycles to let the filter converge
    for(int i = 0; i < 50; ++i)
    {
        filter.predict(Vector<double, 3>::Zero());
        REQUIRE(filter.update(gravity).has_value());
    }

    // After convergence, innovation should be small
    double innov_norm = filter.innovation().norm();
    REQUIRE(innov_norm < 0.1);
    REQUIRE(std::isfinite(innov_norm));
}

TEST_CASE("manifold_ukf create rejects a zero initial quaternion", "[manifold_ukf]")
{
    // Before the fallible factory existed, the constructor normalized the zero
    // quaternion directly and silently produced an all-NaN filter state.
    manifold_ukf_config<double, 3> cfg;
    cfg.q0 = Eigen::Quaternion<double>{0.0, 0.0, 0.0, 0.0};

    auto filter = MukfType::create(simple_rotation_dynamics{}, gravity_meas{}, cfg);
    REQUIRE_FALSE(filter.has_value());
    REQUIRE(filter.error() == filter_error::degenerate_quaternion);
}

TEST_CASE("manifold_ukf create accepts a non-unit nonzero quaternion", "[manifold_ukf]")
{
    // The guard rejects only degeneracy: any finite nonzero quaternion is
    // normalized onto the unit sphere at construction.
    manifold_ukf_config<double, 3> cfg;
    cfg.q0 = Eigen::Quaternion<double>{2.0, 0.0, 0.0, 0.0};

    auto filter = MukfType::create(simple_rotation_dynamics{}, gravity_meas{}, cfg);
    REQUIRE(filter.has_value());
    REQUIRE((filter->attitude().coeffs().array() == cfg.q0.normalized().coeffs().array()).all());
}
