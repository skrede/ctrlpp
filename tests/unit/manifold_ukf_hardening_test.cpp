// What the oracles in this file decide.
//
// The filter carries an attitude quaternion and a three-dimensional tangent
// covariance, and every case feeds it a gravity measurement, which is blind to
// rotation about gravity. The claims here are exact wherever the arithmetic
// allows and derived from the configuration everywhere else:
//
//  * A rejected measurement leaves the attitude, the covariance and the
//    quaternion norm BITWISE unchanged, names its enumerator, and does not
//    degrade the health status.
//  * "Normalizes a non-unit initial quaternion" is asserted as what that phrase
//    MEANS: a filter configured with a non-unit quaternion is BITWISE
//    indistinguishable, in state and in covariance, from one configured with
//    its normalized form, over a run driven by a non-trivial rate and a
//    non-trivial measurement. Finiteness, or even a unit norm, would pass for a
//    filter that normalized to the wrong rotation.
//  * Positive definiteness carries NO negative floor and symmetry is asserted
//    EXACTLY, on the same argument as the multiplicative filter next door.
//  * Convergence is asserted as an exact fixed point at the truth, a no-floor
//    ladder across four checkpoints an order apart, and agreement between an
//    initial tilt about x and one about y -- the rotational symmetry of a
//    gravity measurement. That last one is NOT bitwise here, unlike its
//    multiplicative counterpart, because the sigma-point recombination sums
//    seven terms weighted at about a million, so the two coordinate orders
//    round differently; the disagreement is bounded by that amplification times
//    a counted operation count times the epsilon.
//  * The near-half-turn case IS a convergence case here, and it converges: the
//    tilt falls STRICTLY at every one of two hundred steps from within five
//    hundredths of a radian of a half turn, and keeps falling afterwards.
//  * The hemisphere case asserts the one thing the hemisphere handling changes:
//    a sigma point placed MORE than a half turn from the mean is measured along
//    the SHORT geodesic when the predicted covariance is formed. The whole
//    predicted diagonal is computed here from the strategy's own spread and the
//    configured covariance, and it is the short distance squared over the
//    tangent dimension, not the long one.
//  * The posterior covariance is checked independently on an affine
//    tangent-state measurement at the identity. Diagonal prior, measurement,
//    and noise matrices reduce the Kalman covariance equation to three scalar
//    closed forms, so the oracle does not repeat the filter's matrix update.
//
// What they deliberately do not decide. No convergence RATE is asserted: the
// rate is the filter's own Riccati-like recursion and stating it independently
// would need a second implementation of the estimator under test.
//
// A finding, recorded here because it governs how the hemisphere case is
// written. The filter carries TWO redundant hemisphere mechanisms: an explicit
// antipodal test at three call sites, and the canonicalization inside the
// logarithm those sites feed. Either alone selects the short representative, so
// removing either one alone changes NO result -- measured, with both the whole
// suite still green. Only removing both moves the predicted covariance, from
// 3.0862 to 3.4999 on this fixture, and that is the mutation the case below is
// written to fail. A case asserting only that a sigma point straddles the
// boundary, or only that the mean is unmoved, cannot fail under any of the
// three mutations, which is why neither is the load-bearing assertion here.

#include "hardening_helpers.h"
#include "ctrlpp/estimation/manifold_ukf.h"
#include "ctrlpp/lie/so3.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <utility>

namespace {

struct simple_rotation_dynamics
{
    double dt = 0.01;

    auto operator()(const Eigen::Quaternion<double>& q,
                    const ctrlpp::Vector<double, 3>& omega) const -> Eigen::Quaternion<double>
    {
        ctrlpp::Vector<double, 3> phi = (omega * dt).eval();
        return (q * ctrlpp::so3::exp(phi)).normalized();
    }
};

struct gravity_meas
{
    auto operator()(const Eigen::Quaternion<double>& q) const -> ctrlpp::Vector<double, 3>
    {
        return q.toRotationMatrix().transpose().col(2);
    }
};

struct extreme_measurement
{
    auto operator()(const Eigen::Quaternion<double>&) const
        -> ctrlpp::Vector<double, 3>
    {
        return ctrlpp::Vector<double, 3>::Constant(
            -std::numeric_limits<double>::max());
    }
};

struct diagonal_tangent_measurement
{
    auto operator()(const Eigen::Quaternion<double>& q) const
        -> ctrlpp::Vector<double, 3>
    {
        const auto tangent = ctrlpp::so3::log(q);
        ctrlpp::Vector<double, 3> measured;
        measured << tangent(0), 2.0 * tangent(1), -0.5 * tangent(2);
        return measured;
    }
};

using MukfType = ctrlpp::manifold_ukf<double, 3, simple_rotation_dynamics, gravity_meas>;

// create() is the only construction path and it is fallible, so every
// valid-input site goes through it and asserts success here.
auto make_filter(const ctrlpp::manifold_ukf_config<double, 3>& cfg) -> MukfType
{
    auto created = MukfType::create(simple_rotation_dynamics{}, gravity_meas{}, cfg);
    REQUIRE(created.has_value());
    return *std::move(created);
}

// The rotation angle, read off the quaternion's vector part rather than through
// an arc cosine of a dot product, which saturates to exactly zero below about
// 1e-8 and would hide the tail of every convergence run in this file.
auto rotation_angle(const Eigen::Quaterniond& q) -> double
{
    return 2.0 * std::atan2(q.vec().norm(), std::abs(q.w()));
}

auto gravity_of_identity() -> ctrlpp::Vector<double, 3>
{
    ctrlpp::Vector<double, 3> z;
    z << 0.0, 0.0, 1.0;
    return z;
}

// A step that moves the attitude ends in a renormalization: one correctly
// rounded division of four coefficients by a norm. It does not accumulate,
// because every step renormalizes from scratch. Measured worst deviation over
// the runs below: half of one such rounding.
constexpr int normalize_rounding_ops = 1;
const double unit_norm_budget = normalize_rounding_ops * std::numeric_limits<double>::epsilon();

// The scaled unscented transform recombines seven sigma points with weights
// whose magnitudes are |lambda| / (n + lambda) and 1 / (2 (n + lambda)). At the
// default spread that is about a million, so a quantity of order one is formed
// as a sum of terms of order a million and the cancellation amplifies each
// rounding by that factor. The count is the recombinations one cycle performs:
// the predicted mean's iteration, the predicted covariance, the measurement
// mean, the innovation covariance and the cross-covariance, at five terms each.
constexpr int transform_rounding_ops = 25;

auto sigma_cancellation_amplification() -> double
{
    // lambda = alpha^2 (n + kappa) - n at the strategy's own defaults.
    const ctrlpp::merwe_options<double> opts{};
    const double n = 3.0;
    const double lambda = opts.alpha * opts.alpha * (n + opts.kappa) - n;
    return std::abs(lambda) / (n + lambda);
}

}

TEST_CASE("Manifold UKF non-unit quaternion input normalizes",
          "[manifold_ukf][hardening][negative]")
{
    // "Normalizes" means the filter behaves as though it had been handed the
    // normalized quaternion. A norm assertion alone would pass for a filter
    // that scaled onto the sphere along the wrong direction, and finiteness
    // would pass for one that did nothing at all.
    auto configured = [](const Eigen::Quaterniond& q0) {
        ctrlpp::manifold_ukf_config<double, 3> cfg;
        cfg.Q *= 1e-6;
        cfg.R *= 0.01;
        cfg.q0 = q0;
        return cfg;
    };

    // Doubling every coefficient of the identity: the normalization is a
    // division by exactly two, which is exact in a binary radix, so the seeded
    // attitude is bitwise the identity.
    auto doubled = make_filter(configured(Eigen::Quaterniond(2.0, 0.0, 0.0, 0.0)));
    CHECK(doubled.attitude().coeffs() == Eigen::Quaterniond::Identity().coeffs());

    // And a quaternion whose normalization is not exact, driven by a rate and a
    // measurement that both break every symmetry of the fixture, so an
    // agreement here cannot come from both filters sitting still.
    const Eigen::Quaterniond skew(3.0, 1.0, -2.0, 0.5);
    auto raw = make_filter(configured(skew));
    auto pre_normalized = make_filter(configured(skew.normalized()));

    CHECK(raw.state() == pre_normalized.state());

    ctrlpp::Vector<double, 3> omega;
    omega << 0.05, -0.02, 0.01;
    ctrlpp::Vector<double, 3> z;
    z << 0.1, 0.0, 0.99;

    for(int k = 0; k < 5; ++k)
    {
        raw.predict(omega);
        pre_normalized.predict(omega);
        REQUIRE(raw.update(z).has_value());
        REQUIRE(pre_normalized.update(z).has_value());
        REQUIRE(raw.state() == pre_normalized.state());
        REQUIRE(raw.covariance() == pre_normalized.covariance());
    }

    CHECK(std::abs(raw.attitude().norm() - 1.0) <= unit_norm_budget);
}

TEST_CASE("Manifold UKF rejects correction overflow without committing it",
          "[manifold_ukf][hardening][negative]")
{
    ctrlpp::manifold_ukf_config<double, 3> config;
    using filter_type =
        ctrlpp::manifold_ukf<double, 3,
                            simple_rotation_dynamics,
                            extreme_measurement>;
    auto filter = ctrlpp::test::constructed(filter_type::create(
        simple_rotation_dynamics{}, extreme_measurement{}, config));
    auto const state_before = filter.state();
    auto const covariance_before = filter.covariance();
    auto const innovation_before = filter.innovation();
    auto const attitude_before = filter.attitude();
    auto const health_before = filter.health();

    auto measurement = ctrlpp::Vector<double, 3>::Constant(
        std::numeric_limits<double>::max());
    auto result = filter.update(measurement);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error()
          == ctrlpp::manifold_ukf_update_error::non_finite_result);
    CHECK(filter.state() == state_before);
    CHECK(filter.covariance() == covariance_before);
    CHECK(filter.innovation() == innovation_before);
    CHECK(filter.attitude().coeffs() == attitude_before.coeffs());
    CHECK(filter.health() == health_before);
}

TEST_CASE("Manifold UKF NaN rotation measurement is rejected without touching the estimate",
          "[manifold_ukf][hardening][negative]")
{
    ctrlpp::manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-6;
    cfg.R *= 0.01;

    auto filter = make_filter(cfg);
    // Stepped only with the valid measurement, never with the poisoned one, so
    // it says what the filter would have carried had the bad sample never
    // arrived.
    auto reference = make_filter(cfg);

    ctrlpp::Vector<double, 3> omega = ctrlpp::Vector<double, 3>::Zero();
    filter.predict(omega);
    reference.predict(omega);

    // Snapshot immediately before the poisoned step.
    const ctrlpp::Vector<double, 4> state_before = filter.state();
    const ctrlpp::Matrix<double, 3, 3> P_before = filter.covariance();

    ctrlpp::Vector<double, 3> z_bad;
    z_bad << std::numeric_limits<double>::quiet_NaN(), 0.0, 1.0;

    const auto rejected = filter.update(z_bad);

    REQUIRE_FALSE(rejected.has_value());
    REQUIRE(rejected.error() == ctrlpp::manifold_ukf_update_error::non_finite_measurement);

    // Exact comparison, not a tolerance: a rejected step performs no arithmetic
    // on the carried estimate at all, so bitwise equality is the contract and a
    // tolerance would admit a step that partially ran.
    CHECK(filter.state() == state_before);
    // The covariance was already measurement-independent before the guard
    // existed -- the tangent reduction P - K*S*K^T and its reset Jacobian are
    // built from the sigma points, the gain and the correction, never from z --
    // so this half of the invariant is structural. The state half is what the
    // guard adds.
    CHECK(filter.covariance() == P_before);
    // On a manifold the attitude quaternion's unit norm is an extra invariant,
    // and the rejection preserves it exactly rather than approximately.
    CHECK(filter.attitude().norm() == 1.0);
    // A rejection describes the sample, not the filter: nothing was mutated, so
    // the filter is not degraded and must not report that it is.
    CHECK(filter.health() == ctrlpp::manifold_ukf_health::ok);

    // The poison did not latch: the next valid step produces exactly what it
    // would have produced had the poisoned step never been attempted.
    ctrlpp::Vector<double, 3> z_good;
    z_good << 0.0, 0.0, 1.0;
    REQUIRE(filter.update(z_good).has_value());
    REQUIRE(reference.update(z_good).has_value());

    CHECK(filter.state() == reference.state());
    CHECK(filter.covariance() == reference.covariance());
}

TEST_CASE("Manifold UKF covariance stays symmetric positive definite over 1000 steps",
          "[manifold_ukf][hardening][stability]")
{
    ctrlpp::manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-6;
    cfg.R *= 0.01;

    auto filter = make_filter(cfg);

    ctrlpp::Vector<double, 3> omega = ctrlpp::Vector<double, 3>::Zero();

    for(int k = 0; k < 1000; ++k)
    {
        filter.predict(omega);
        REQUIRE(filter.update(gravity_of_identity()).has_value());

        // Symmetry is EXACT: both the predicted covariance and the posterior
        // reduction pass through the symmetrizing helper, so the two triangles
        // hold the same bits.
        REQUIRE(filter.covariance() == filter.covariance().transpose());

        // Positive definiteness with NO negative floor. A floor admits an
        // indefinite covariance in a case named for definiteness.
        Eigen::SelfAdjointEigenSolver<ctrlpp::Matrix<double, 3, 3>> eigsolver(filter.covariance());
        REQUIRE(eigsolver.eigenvalues().minCoeff() > 0.0);
    }

    // Every geodesic mean converged within its budget, which is what makes the
    // predicted attitudes above means rather than last iterates.
    CHECK(filter.health() == ctrlpp::manifold_ukf_health::ok);
}

TEST_CASE("Manifold UKF posterior covariance matches an affine tangent-state closed form",
          "[manifold_ukf][hardening][covariance]")
{
    using filter_type =
        ctrlpp::manifold_ukf<double, 3, simple_rotation_dynamics,
                            diagonal_tangent_measurement>;

    constexpr std::array<double, 3> prior_variances{0.04, 0.09, 0.16};
    constexpr std::array<double, 3> measurement_gains{1.0, 2.0, -0.5};
    constexpr std::array<double, 3> measurement_variances{0.01, 0.04, 0.25};

    ctrlpp::manifold_ukf_config<double, 3> cfg;
    cfg.P0 = ctrlpp::Matrix<double, 3, 3>::Zero();
    cfg.R = ctrlpp::Matrix<double, 3, 3>::Zero();
    for(int i = 0; i < 3; ++i)
    {
        cfg.P0(i, i) = prior_variances[static_cast<std::size_t>(i)];
        cfg.R(i, i) = measurement_variances[static_cast<std::size_t>(i)];
    }

    // alpha = 1 and kappa = 0 make lambda zero. The center sigma point then
    // contributes no covariance, and each opposite pair carries one sixth of
    // the tangent spread. That makes an affine tangent measurement reproduce
    // P, H P H^T, and P H^T without the million-fold cancellation of the
    // default spread.
    const ctrlpp::merwe_options<double> wide{
        .alpha = 1.0, .beta = 0.0, .kappa = 0.0};
    auto strategy = ctrlpp::so3_merwe_sigma_points<double>::try_create(wide);
    REQUIRE(strategy.has_value());

    auto filter = ctrlpp::test::constructed(filter_type::create(
        simple_rotation_dynamics{}, diagonal_tangent_measurement{}, cfg,
        *strategy));

    const auto zero_measurement = ctrlpp::Vector<double, 3>::Zero();
    REQUIRE(filter.update(zero_measurement).has_value());

    // Opposite affine measurements cancel exactly on this fixture, so the
    // correction and reset Jacobian are the identity. This premise is asserted:
    // if it moves, the scalar posterior below is no longer the applicable
    // closed form and must not be allowed to pass accidentally.
    REQUIRE(filter.innovation() == zero_measurement);
    REQUIRE(filter.attitude().coeffs()
            == Eigen::Quaterniond::Identity().coeffs());

    // A conservative forward-error budget counted from the six non-central
    // sigma points: each exp/log tangent round trip and diagonal measurement
    // uses at most 27 rounded operations (162); their weighted measurement and
    // cross-covariance accumulations use at most 42 more; and the diagonal
    // innovation solve, covariance reduction, scalar oracle, and comparison use
    // at most 52. Every contribution is bounded at the largest covariance
    // operand's scale, for 256 operations in total.
    constexpr int posterior_covariance_rounding_ops = 256;
    const double eps = std::numeric_limits<double>::epsilon();

    for(int i = 0; i < 3; ++i)
    {
        const auto index = static_cast<std::size_t>(i);
        const double prior = prior_variances[index];
        const double gain = measurement_gains[index];
        const double noise = measurement_variances[index];

        // Scalar Kalman posterior for z_i = gain * x_i:
        //
        //   p+ = p - p^2 gain^2 / (gain^2 p + r)
        //      = p r / (gain^2 p + r).
        //
        // This contains neither the filter's K*S*K' expression nor a matrix
        // factorization, so a multiplicative error in that reduction cannot
        // reproduce the oracle.
        const double expected = prior * noise
                                / (gain * gain * prior + noise);
        const double scale = std::max({prior, noise, expected});
        const double tolerance =
            static_cast<double>(posterior_covariance_rounding_ops) * eps
            * scale;

        CAPTURE(i, prior, gain, noise, expected,
                filter.covariance()(i, i), tolerance);
        CHECK(std::abs(filter.covariance()(i, i) - expected) <= tolerance);
    }

    const double off_diagonal_scale = prior_variances.back();
    const double off_diagonal_tolerance =
        static_cast<double>(posterior_covariance_rounding_ops) * eps
        * off_diagonal_scale;
    CHECK(std::abs(filter.covariance()(0, 1)) <= off_diagonal_tolerance);
    CHECK(std::abs(filter.covariance()(0, 2)) <= off_diagonal_tolerance);
    CHECK(std::abs(filter.covariance()(1, 2)) <= off_diagonal_tolerance);
    CHECK(filter.covariance() == filter.covariance().transpose());
}

TEST_CASE("Manifold UKF attitude converges for slow rotation",
          "[manifold_ukf][hardening][convergence]")
{
    constexpr double initial_tilt = 0.3;

    auto configured = [](const Eigen::Quaterniond& q0) {
        ctrlpp::manifold_ukf_config<double, 3> cfg;
        cfg.Q *= 1e-6;
        cfg.R *= 0.01;
        cfg.q0 = q0;
        return cfg;
    };

    // The truth is an EXACT fixed point: seeded at the identity and fed the
    // gravity of the identity, the filter must not move at all.
    {
        auto at_truth = make_filter(configured(Eigen::Quaterniond::Identity()));
        const ctrlpp::Vector<double, 4> seeded = at_truth.state();
        for(int k = 0; k < 500; ++k)
        {
            at_truth.predict(ctrlpp::Vector<double, 3>::Zero());
            REQUIRE(at_truth.update(gravity_of_identity()).has_value());
            REQUIRE(at_truth.state() == seeded);
        }
    }

    // The tilt keeps falling across checkpoints an order apart. A filter that
    // converges two thirds of the way and stalls passes any single threshold
    // below its own starting error and fails this.
    auto filter = make_filter(configured(
        Eigen::Quaterniond(Eigen::AngleAxisd(initial_tilt, Eigen::Vector3d::UnitY()))));

    double at_100 = 0.0;
    double at_500 = 0.0;
    double at_1000 = 0.0;
    double at_2000 = 0.0;

    for(int k = 1; k <= 2000; ++k)
    {
        filter.predict(ctrlpp::Vector<double, 3>::Zero());
        REQUIRE(filter.update(gravity_of_identity()).has_value());
        REQUIRE(std::abs(filter.attitude().norm() - 1.0) <= unit_norm_budget);

        const double tilt = rotation_angle(filter.attitude());
        if(k == 100)
            at_100 = tilt;
        if(k == 500)
            at_500 = tilt;
        if(k == 1000)
            at_1000 = tilt;
        if(k == 2000)
            at_2000 = tilt;
    }

    CHECK(at_100 < initial_tilt);
    CHECK(at_500 < at_100);
    CHECK(at_1000 < at_500);
    CHECK(at_2000 < at_1000);

    // A gravity measurement is invariant under rotation about gravity, so an
    // initial tilt about x and one about y are the same problem in different
    // coordinates. The two trajectories agree to the level the sigma-point
    // recombination allows: a quantity of order one assembled from terms of
    // order a million, so each rounding is amplified by that ratio.
    const double symmetry_budget = transform_rounding_ops * sigma_cancellation_amplification()
        * std::numeric_limits<double>::epsilon();

    auto about_x = make_filter(configured(
        Eigen::Quaterniond(Eigen::AngleAxisd(initial_tilt, Eigen::Vector3d::UnitX()))));
    auto about_y = make_filter(configured(
        Eigen::Quaterniond(Eigen::AngleAxisd(initial_tilt, Eigen::Vector3d::UnitY()))));

    for(int k = 0; k < 500; ++k)
    {
        about_x.predict(ctrlpp::Vector<double, 3>::Zero());
        about_y.predict(ctrlpp::Vector<double, 3>::Zero());
        REQUIRE(about_x.update(gravity_of_identity()).has_value());
        REQUIRE(about_y.update(gravity_of_identity()).has_value());

        const double tilt_x = rotation_angle(about_x.attitude());
        const double tilt_y = rotation_angle(about_y.attitude());
        REQUIRE(std::abs(tilt_x - tilt_y) <= symmetry_budget * std::max(tilt_x, tilt_y));
    }
}

TEST_CASE("Manifold UKF converges through a near-half-turn tilt",
          "[manifold_ukf][hardening][robustness]")
{
    // The fixture starts within five hundredths of a radian of a half turn
    // ABOUT AN AXIS PERPENDICULAR TO GRAVITY, so the whole error is observable
    // and the question the fixture poses -- does the filter come back through
    // the antipodal region -- has an answer. It does, and monotonically.
    const Eigen::Quaterniond q_init(Eigen::AngleAxisd(
        std::numbers::pi - 0.05, Eigen::Vector3d::UnitX()));

    ctrlpp::manifold_ukf_config<double, 3> cfg;
    cfg.Q *= 1e-4;
    cfg.R *= 0.1;
    cfg.q0 = q_init;

    auto filter = make_filter(cfg);

    const double initial_tilt = rotation_angle(filter.attitude());
    double previous = initial_tilt;

    for(int k = 0; k < 200; ++k)
    {
        filter.predict(ctrlpp::Vector<double, 3>::Zero());
        REQUIRE(filter.update(gravity_of_identity()).has_value());
        REQUIRE(std::abs(filter.attitude().norm() - 1.0) <= unit_norm_budget);

        // Strictly, at every step, from a tilt of more than three radians. No
        // threshold is chosen: a filter that stalls, oscillates or diverges
        // anywhere in the antipodal region fails at the step where it does.
        const double tilt = rotation_angle(filter.attitude());
        REQUIRE(tilt < previous);
        previous = tilt;
    }

    const double after_200 = previous;
    CHECK(after_200 < initial_tilt);

    // And it does not stall afterwards.
    for(int k = 0; k < 200; ++k)
    {
        filter.predict(ctrlpp::Vector<double, 3>::Zero());
        REQUIRE(filter.update(gravity_of_identity()).has_value());
    }
    CHECK(rotation_angle(filter.attitude()) < after_200);
    CHECK(filter.health() == ctrlpp::manifold_ukf_health::ok);
}

TEST_CASE("Manifold UKF measures sigma points past the half turn along the shorter geodesic",
          "[manifold_ukf][hardening][coverage]")
{
    // The property the hemisphere handling delivers, and the only one that is
    // observable: a sigma point placed MORE than a half turn from the mean is
    // the same rotation as one placed the short way round in the opposite
    // sense, and the predicted covariance -- which is a spread measured on the
    // manifold -- must be built from the SHORT distance. A filter that measured
    // the long way would over-report its own uncertainty.
    //
    // The fixture has to be built for this. At the strategy's default spread
    // the tangent offsets are the spread parameter times the square root of the
    // covariance, and the spread parameter is a thousandth, so no covariance in
    // a sane range places a sigma point past a half turn. The widened spread
    // below does, and that is asserted rather than assumed.
    const ctrlpp::merwe_options<double> wide{.alpha = 1.0, .beta = 2.0, .kappa = 0.0};
    auto strategy_result = ctrlpp::so3_merwe_sigma_points<double>::try_create(wide);
    REQUIRE(strategy_result.has_value());
    const auto strategy = *strategy_result;

    constexpr double tangent_dimension = 3.0;
    constexpr double initial_variance = 3.5;
    constexpr double process_variance = 1e-6;

    const Eigen::Quaterniond q_near_pi(Eigen::AngleAxisd(3.0, Eigen::Vector3d::UnitZ()));
    const Eigen::Quaterniond antipode(-q_near_pi.w(), -q_near_pi.x(), -q_near_pi.y(), -q_near_pi.z());

    auto configured = [&](const Eigen::Quaterniond& q0) {
        ctrlpp::manifold_ukf_config<double, 3> cfg;
        cfg.q0 = q0;
        cfg.P0 = ctrlpp::Matrix<double, 3, 3>::Identity() * initial_variance;
        cfg.Q = ctrlpp::Matrix<double, 3, 3>::Identity() * process_variance;
        cfg.R *= 0.1;
        return cfg;
    };

    auto direct = MukfType::create(simple_rotation_dynamics{}, gravity_meas{},
                                   configured(q_near_pi), strategy);
    auto flipped = MukfType::create(simple_rotation_dynamics{}, gravity_meas{},
                                    configured(antipode), strategy);
    REQUIRE(direct.has_value());
    REQUIRE(flipped.has_value());

    // The premise: the generated set really does straddle the boundary. A
    // negative scalar product with the mean is exactly the condition the
    // hemisphere handling exists for, and six of the seven points meet it here.
    {
        const auto generated = strategy.generate(direct->attitude(), direct->covariance());
        int straddling = 0;
        for(const auto& point : generated.points)
        {
            if(point.dot(direct->attitude()) < 0.0)
                ++straddling;
        }
        REQUIRE(straddling > 0);
    }

    // With no commanded rate the propagation is the identity, so the predicted
    // covariance is the sigma set's own spread about its centre plus the
    // process noise, and every factor of it is fixed by the configuration:
    //
    //   offset radius   = sqrt(n + lambda) * sqrt(P0)   = sqrt(3 * 3.5)
    //                   = 3.2404 radians, which EXCEEDS a half turn
    //   short distance  = 2 pi - offset radius          = 3.0428 radians
    //   diagonal weight = 2 * Wc_i = 1 / (n + lambda)   = 1/3
    //
    // so each diagonal entry is the square of the short distance over three,
    // plus the process variance. The two offsets of a given axis contribute
    // only to that axis's diagonal, the centre point contributes a zero tangent
    // vector, and lambda is exactly zero at this spread so the centre weight
    // drops out of the mean.
    const double offset_radius = std::sqrt(tangent_dimension * initial_variance);
    REQUIRE(offset_radius > std::numbers::pi);
    const double short_distance = 2.0 * std::numbers::pi - offset_radius;
    const double predicted_variance = short_distance * short_distance / tangent_dimension + process_variance;

    // Nine roundings: the square root forming the offset radius (1), the
    // Cholesky factor's own square root (2), the scaling of the tangent basis
    // (3), the exponential's half angle, cosine and cardinal sine (4-6), the
    // logarithm's arc tangent and its scaling (7-8), and the weighted
    // accumulation of the outer product (9).
    constexpr int predicted_covariance_ops = 9;
    const double budget = predicted_covariance_ops * std::numeric_limits<double>::epsilon() * predicted_variance;

    {
        auto centred = *MukfType::create(simple_rotation_dynamics{}, gravity_meas{},
                                         configured(q_near_pi), strategy);
        const Eigen::Quaterniond before = centred.attitude();
        centred.predict(ctrlpp::Vector<double, 3>::Zero());

        // The mean is unmoved: each antipodal pair cancels whichever
        // representative is chosen, so this half is true with or without the
        // hemisphere handling and is asserted as a precondition, not as the
        // property under test.
        CHECK(centred.attitude().coeffs() == before.coeffs());

        // This is the property under test. Measured with both hemisphere
        // mechanisms removed -- the explicit antipodal test in the covariance
        // sum and the canonicalization inside the logarithm -- each diagonal
        // becomes the LONG distance squared over three, 3.4999 against 3.0862,
        // and this assertion fails. Removing either one alone changes nothing,
        // because the other still selects the short representative.
        for(int i = 0; i < 3; ++i)
            CHECK(std::abs(centred.covariance()(i, i) - predicted_variance) <= budget);
    }

    // A quaternion and its negation are the same rotation, so the two filters
    // must agree. Recorded honestly: this holds by the conjugation structure of
    // the recursion -- every tangent quantity is formed from q_mean^{-1} q_i,
    // where a sign on both factors cancels -- and it therefore survives the
    // removal of the hemisphere handling. It is kept because it pins a real
    // equivariance the filter must not lose, not because it covers the branch.
    Eigen::Vector3d omega(0.1, 0.0, 0.5);
    const Eigen::Vector3d z = q_near_pi.toRotationMatrix().transpose().col(2);

    for(int k = 0; k < 50; ++k)
    {
        direct->predict(omega);
        flipped->predict(omega);
        REQUIRE(direct->update(z).has_value());
        REQUIRE(flipped->update(z).has_value());

        REQUIRE(direct->attitude().coeffs() == -flipped->attitude().coeffs());
        REQUIRE(direct->covariance() == flipped->covariance());
    }
}
