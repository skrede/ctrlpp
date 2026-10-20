// What the oracles in this file decide.
//
// The linear filter is the one estimator here whose answer is fixed in closed
// form at every operating point the cases use, so the file asserts closed forms
// and typed rejections, never finiteness:
//
//  * The typed rejections: each non-finite, structurally asymmetric or
//    indefinite configuration field by name, and a non-finite measurement, the
//    latter leaving state and covariance BITWISE unchanged.
//  * The deadbeat identity when the measurement noise is zero -- the corrected
//    estimate IS the measurement in the observed coordinate and the posterior
//    covariance there is exactly zero, both bitwise.
//  * The covariance contraction with zero process noise, asserted as the exact
//    determinant identity det(P_post) * S = det(P_pred) * R, which holds because
//    the state matrix here has unit determinant and the optimal gain contracts
//    the covariance volume by exactly R/S. This is the theorem; a monotone
//    trace is NOT, and although the trace is observed to decrease monotonically
//    for this data it is not asserted, because a congruence by a state matrix
//    with a Jordan block can raise it for other data.
//  * Estimation error propagated exactly. The measurement streams below are
//    noiseless and exactly consistent with the model, so the error obeys the
//    homogeneous recursion e(k+1) = (I - K(k) C) A e(k) with the filter's own
//    time-varying gain. The oracle is that propagated error, recomputed in the
//    test from the filter's own reported covariance -- not a fitted final
//    threshold and not a comparison against the filter's reported sigma.
//  * The unobservable mode's covariance growth in closed form, and its estimate
//    fixed point in closed form.
//
// What they deliberately do not decide. Nothing here asserts that the estimate
// converges to the truth in a statistical sense, because no case here injects
// measurement noise; the streams are deterministic and the error is a transient
// with a closed form, which is a stronger thing to assert and a different one.
// Nothing here asserts optimality of the gain against an independent optimizer
// either -- the determinant identity and the propagated error both read the
// gain the filter chose.

#include "hardening_helpers.h"
#include "ctrlpp/estimation/kalman.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

constexpr double kf_eps = std::numeric_limits<double>::epsilon();

// Rounded operations along the longest chain through one filter step for a
// two-state, one-output problem. Enumerated rather than chosen:
//   state propagation      x = A x + B u          5  (two products and a sum per
//                                                     component, plus the input
//                                                     product and its sum)
//   covariance propagation P = A P A' + Q         7  (two contractions of length
//                                                     two, three operations each,
//                                                     plus the process-noise sum)
//   innovation             y = z - C x            4
//   innovation covariance  S = C P C' + R         7
//   gain                   K from a QR solve      5  (the contraction C P, three
//                                                     operations, plus a backward
//                                                     error bounded by twice the
//                                                     output dimension)
//   state correction       x = x + K y            2
//   Joseph covariance      (I-KC)P(I-KC)' + K R K' 13 (the factor, two
//                                                     contractions, the noise
//                                                     term, the sum and the
//                                                     symmetrization)
// Every operation is counted whether or not it actually rounds, so the total
// bounds the accumulated error from above rather than describing it tightly --
// which is what a budget requires.
constexpr int kf_step_ops = 43;

// Rounded operations reaching the determinant identity below. The Joseph chain
// above contributes thirteen and the gain solve five; each two-by-two
// determinant costs two products and one difference, and the innovation
// covariance costs three more -- twenty-four in all. The scale is NOT the
// determinant: a two-by-two determinant is a difference of two products that
// can cancel almost exactly, so the rounding is set by the magnitudes of those
// products, which the test forms explicitly.
constexpr int kf_determinant_ops = 24;

// The symmetric eigensolver's backward perturbation is bounded by a small
// multiple of the number of matrix entries it touches, times epsilon, times the
// norm of the matrix -- the Weyl perturbation bound. For a two-by-two that is
// four entries.
constexpr int kf_eig_ops = 4;

auto make_const_velocity_system()
{
    ctrlpp::discrete_state_space<double, 2, 1, 1> sys;
    double dt = 0.1;
    sys.A << 1.0, dt, 0.0, 1.0;
    sys.B << 0.5 * dt * dt, dt;
    sys.C << 1.0, 0.0;
    sys.D << 0.0;
    return sys;
}

/// @brief The gain the filter must have chosen at its current predicted
/// covariance, and the closed-loop error map that gain produces.
///
/// Recomputed in the test from the filter's OWN reported predicted covariance,
/// so it reads what the filter did rather than assuming what it should have
/// done. The estimation error of a noiseless, model-consistent measurement
/// stream obeys e(k+1) = (I - K C) A e(k) exactly, which is what makes this the
/// oracle for the convergence cases.
struct step_gain
{
    Eigen::Vector2d K;
    Eigen::Matrix2d closed_loop;
    double innovation_covariance;
};

auto gain_at(const ctrlpp::discrete_state_space<double, 2, 1, 1>& sys,
             const Eigen::Matrix2d& P_pred, double R) -> step_gain
{
    const double S = (sys.C * P_pred * sys.C.transpose())(0, 0) + R;
    const Eigen::Vector2d K = (P_pred * sys.C.transpose()).eval() / S;
    return {K, ((Eigen::Matrix2d::Identity() - K * sys.C) * sys.A).eval(), S};
}

}

TEST_CASE("Kalman rejects each non-finite configuration field by name", "[kalman][hardening][negative]")
{
    auto sys = make_const_velocity_system();

    // The failure this prevents: an infinite Q makes the first predict give
    // P = A P A' + Inf = Inf, the gain solve is posed against S = C P C' + R =
    // Inf and yields Inf/Inf = NaN, and the corrected state follows. The caller
    // would see a non-finite estimate from a filter it configured, with nothing
    // naming the field that was wrong.
    SECTION("process noise")
    {
        ctrlpp::kalman_config<double, 2, 1, 1> cfg{};
        cfg.Q = ctrlpp::test::inf_matrix<double, 2, 2>();
        const auto rejected = ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::non_finite_process_noise);
    }

    SECTION("measurement noise")
    {
        ctrlpp::kalman_config<double, 2, 1, 1> cfg{};
        cfg.R = ctrlpp::test::nan_matrix<double, 1, 1>();
        const auto rejected = ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::non_finite_measurement_noise);
    }

    SECTION("initial state")
    {
        ctrlpp::kalman_config<double, 2, 1, 1> cfg{};
        cfg.x0 = ctrlpp::test::nan_vector<double, 2>();
        const auto rejected = ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::non_finite_initial_state);
    }

    SECTION("initial covariance")
    {
        ctrlpp::kalman_config<double, 2, 1, 1> cfg{};
        cfg.P0 = ctrlpp::test::inf_matrix<double, 2, 2>();
        const auto rejected = ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::non_finite_initial_covariance);
    }
}

TEST_CASE("Kalman refuses a configuration matrix that is not a covariance", "[kalman][hardening][negative]")
{
    auto sys = make_const_velocity_system();
    const Eigen::Matrix2d indefinite = Eigen::Vector2d{1.0, -1.0}.asDiagonal();

    SECTION("process noise")
    {
        ctrlpp::kalman_config<double, 2, 1, 1> cfg{};
        cfg.Q = indefinite;
        const auto rejected = ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::indefinite_process_noise);
    }

    SECTION("measurement noise")
    {
        ctrlpp::kalman_config<double, 2, 1, 1> cfg{};
        cfg.R << -1.0;
        const auto rejected = ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::indefinite_measurement_noise);
    }

    SECTION("initial covariance")
    {
        ctrlpp::kalman_config<double, 2, 1, 1> cfg{};
        cfg.P0 = indefinite;
        const auto rejected = ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::indefinite_initial_covariance);
    }

    SECTION("an earlier field's finiteness is reported ahead of a later field's definiteness")
    {
        ctrlpp::kalman_config<double, 2, 1, 1> cfg{};
        cfg.Q = ctrlpp::test::inf_matrix<double, 2, 2>();
        cfg.P0 = indefinite;
        const auto rejected = ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::non_finite_process_noise);
    }
}

TEST_CASE("Kalman refuses a structurally asymmetric configuration matrix by name", "[kalman][hardening][negative]")
{
    auto sys = make_const_velocity_system();
    Eigen::Matrix2d asymmetric;
    asymmetric << 1.0, 0.5, 0.4, 1.0;

    SECTION("process noise")
    {
        ctrlpp::kalman_config<double, 2, 1, 1> cfg{};
        cfg.Q = asymmetric;
        const auto rejected = ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::asymmetric_process_noise);
    }

    SECTION("measurement noise of a two-output filter")
    {
        ctrlpp::discrete_state_space<double, 2, 1, 2> two_output;
        two_output.A = sys.A;
        two_output.B = sys.B;
        two_output.C.setIdentity();
        two_output.D.setZero();
        ctrlpp::kalman_config<double, 2, 1, 2> cfg{};
        cfg.R = asymmetric;
        const auto rejected = ctrlpp::kalman_filter<double, 2, 1, 2>::create(two_output, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::asymmetric_measurement_noise);
    }

    SECTION("initial covariance")
    {
        ctrlpp::kalman_config<double, 2, 1, 1> cfg{};
        cfg.P0 = asymmetric;
        const auto rejected = ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::asymmetric_initial_covariance);
    }

    SECTION("a matrix's finiteness is reported ahead of its asymmetry")
    {
        ctrlpp::kalman_config<double, 2, 1, 1> cfg{};
        cfg.Q = asymmetric;
        cfg.Q(1, 1) = std::numeric_limits<double>::infinity();
        const auto rejected = ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::non_finite_process_noise);
    }

    SECTION("an earlier field's asymmetry is reported ahead of a later field's definiteness")
    {
        ctrlpp::kalman_config<double, 2, 1, 1> cfg{};
        cfg.Q = asymmetric;
        cfg.P0 = Eigen::Vector2d{1.0, -1.0}.asDiagonal();
        const auto rejected = ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg);
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::filter_error::asymmetric_process_noise);
    }
}

TEST_CASE("Kalman stores a drift-asymmetric configuration matrix as its midpoint symmetrization", "[kalman][hardening]")
{
    auto sys = make_const_velocity_system();
    ctrlpp::kalman_config<double, 2, 1, 1> drifted{};
    drifted.Q << 0.02, 0.005, std::nextafter(0.005, 1.0), 0.01;
    drifted.P0 << 2.0, 0.3, std::nextafter(0.3, 0.0), 1.0;

    ctrlpp::kalman_config<double, 2, 1, 1> symmetric = drifted;
    symmetric.Q(0, 1) = symmetric.Q(1, 0) = drifted.Q(0, 1) + (drifted.Q(1, 0) - drifted.Q(0, 1)) / 2.0;
    symmetric.P0(0, 1) = symmetric.P0(1, 0) = drifted.P0(0, 1) + (drifted.P0(1, 0) - drifted.P0(0, 1)) / 2.0;

    auto filter = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, drifted));
    auto twin = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, symmetric));
    CHECK((filter.covariance().array() == symmetric.P0.array()).all());

    const Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Constant(0.3);
    for(int k = 0; k < 10; ++k)
    {
        const Eigen::Matrix<double, 1, 1> z = Eigen::Matrix<double, 1, 1>::Constant(1.0 + 0.05 * k);
        filter.predict(u);
        twin.predict(u);
        REQUIRE(filter.update(z).has_value());
        REQUIRE(twin.update(z).has_value());
    }
    CHECK((filter.state().array() == twin.state().array()).all());
    CHECK((filter.covariance().array() == twin.covariance().array()).all());
}

// At this scale the relative floor 2 eps max|Q_ij| underflows to zero, so the
// eigenvalue -denorm_min is judged against the absolute floor 2 denorm_min.
TEST_CASE("Kalman accepts a subnormal-scale process noise within the underflow floor", "[kalman][hardening]")
{
    auto sys = make_const_velocity_system();
    const double dm = std::numeric_limits<double>::denorm_min();
    ctrlpp::kalman_config<double, 2, 1, 1> cfg{};
    cfg.Q = Eigen::Vector2d{std::ldexp(dm, 20), -dm}.asDiagonal();
    CHECK(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg).has_value());
}

TEST_CASE("Kalman accepts a singular positive semidefinite configuration", "[kalman][hardening]")
{
    auto sys = make_const_velocity_system();
    ctrlpp::kalman_config<double, 2, 1, 1> cfg{};

    SECTION("rank-one process noise") { cfg.Q << 1.0, 2.0, 2.0, 4.0; }

    SECTION("noiseless measurement channel") { cfg.R << 0.0; }

    SECTION("point-mass prior") { cfg.P0.setZero(); }

    // The floor is 2 eps max|P0_ij| for this 2x2: an eigenvalue at half of it
    // below zero is indistinguishable from zero for the symmetric eigensolver.
    SECTION("an eigenvalue inside the eigensolver's backward error")
    {
        cfg.P0 = Eigen::Vector2d{1.0, -kf_eps}.asDiagonal();
    }

    const auto accepted = ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg);
    REQUIRE(accepted.has_value());
    CHECK(accepted->covariance() == cfg.P0);
    CHECK(accepted->state() == cfg.x0);
}

TEST_CASE("Kalman setters refuse a matrix that is not a covariance and assign nothing", "[kalman][hardening][negative]")
{
    auto sys = make_const_velocity_system();
    ctrlpp::kalman_config<double, 2, 1, 1> cfg{};
    cfg.Q = Eigen::Matrix2d::Identity() * 0.01;
    auto filter = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg));
    auto twin = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg));

    const Eigen::Matrix2d indefinite = Eigen::Vector2d{1.0, -1.0}.asDiagonal();
    const Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();
    const Eigen::Matrix<double, 1, 1> z = Eigen::Matrix<double, 1, 1>::Constant(1.0);

    // Each refused call carries one valid matrix that differs from the
    // filter's, so an assignment of either one would part the twins.
    SECTION("set_noise with an indefinite Q")
    {
        const auto refused = filter.set_noise(indefinite, Eigen::Matrix<double, 1, 1>::Constant(4.0));
        REQUIRE_FALSE(refused.has_value());
        CHECK(refused.error() == ctrlpp::filter_error::indefinite_process_noise);
    }

    SECTION("set_noise with an indefinite R")
    {
        const auto refused = filter.set_noise(Eigen::Matrix2d::Identity() * 2.0, Eigen::Matrix<double, 1, 1>::Constant(-1.0));
        REQUIRE_FALSE(refused.has_value());
        CHECK(refused.error() == ctrlpp::filter_error::indefinite_measurement_noise);
    }

    SECTION("set_noise with a non-finite R")
    {
        const auto refused = filter.set_noise(Eigen::Matrix2d::Identity() * 2.0, ctrlpp::test::nan_matrix<double, 1, 1>());
        REQUIRE_FALSE(refused.has_value());
        CHECK(refused.error() == ctrlpp::filter_error::non_finite_measurement_noise);
    }

    SECTION("reset_covariance with an indefinite P0")
    {
        filter.predict(u);
        twin.predict(u);
        const auto refused = filter.reset_covariance(indefinite);
        REQUIRE_FALSE(refused.has_value());
        CHECK(refused.error() == ctrlpp::filter_error::indefinite_initial_covariance);
        CHECK(filter.covariance() == twin.covariance());
    }

    SECTION("reset_covariance with a non-finite P0")
    {
        const auto refused = filter.reset_covariance(ctrlpp::test::inf_matrix<double, 2, 2>());
        REQUIRE_FALSE(refused.has_value());
        CHECK(refused.error() == ctrlpp::filter_error::non_finite_initial_covariance);
    }

    SECTION("set_noise with a structurally asymmetric Q")
    {
        Eigen::Matrix2d asymmetric;
        asymmetric << 1.0, 0.5, 0.4, 1.0;
        const auto refused = filter.set_noise(asymmetric, Eigen::Matrix<double, 1, 1>::Constant(4.0));
        REQUIRE_FALSE(refused.has_value());
        CHECK(refused.error() == ctrlpp::filter_error::asymmetric_process_noise);
    }

    SECTION("reset_covariance with a structurally asymmetric P0")
    {
        filter.predict(u);
        twin.predict(u);
        Eigen::Matrix2d asymmetric;
        asymmetric << 1.0, 0.5, 0.4, 1.0;
        const auto refused = filter.reset_covariance(asymmetric);
        REQUIRE_FALSE(refused.has_value());
        CHECK(refused.error() == ctrlpp::filter_error::asymmetric_initial_covariance);
        CHECK(filter.covariance() == twin.covariance());
    }

    for(int k = 0; k < 5; ++k)
    {
        filter.predict(u);
        twin.predict(u);
        REQUIRE(filter.update(z).has_value());
        REQUIRE(twin.update(z).has_value());
    }
    CHECK((filter.state().array() == twin.state().array()).all());
    CHECK((filter.covariance().array() == twin.covariance().array()).all());
    CHECK(filter.is_steady_state() == twin.is_steady_state());
}

TEST_CASE("Kalman setters accept a covariance as an assignment", "[kalman][hardening]")
{
    auto sys = make_const_velocity_system();
    ctrlpp::kalman_config<double, 2, 1, 1> cfg{};
    ctrlpp::kalman_config<double, 2, 1, 1> replaced = cfg;
    replaced.Q << 1.0, 2.0, 2.0, 4.0;
    replaced.R << 0.25;
    auto filter = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg));
    auto twin = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, replaced));

    REQUIRE(filter.set_noise(replaced.Q, replaced.R).has_value());

    const Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();
    const Eigen::Matrix<double, 1, 1> z = Eigen::Matrix<double, 1, 1>::Constant(1.0);
    for(int k = 0; k < 5; ++k)
    {
        filter.predict(u);
        twin.predict(u);
        REQUIRE(filter.update(z).has_value());
        REQUIRE(twin.update(z).has_value());
    }
    CHECK((filter.state().array() == twin.state().array()).all());
    CHECK((filter.covariance().array() == twin.covariance().array()).all());

    REQUIRE(filter.reset_covariance(Eigen::Matrix2d::Zero()).has_value());
    CHECK(filter.covariance() == Eigen::Matrix2d::Zero());
}

TEST_CASE("Kalman setters assign a drift-asymmetric matrix as its midpoint symmetrization", "[kalman][hardening]")
{
    auto sys = make_const_velocity_system();
    ctrlpp::kalman_config<double, 2, 1, 1> cfg{};
    ctrlpp::kalman_config<double, 2, 1, 1> replaced = cfg;
    replaced.Q << 0.04, 0.01, std::nextafter(0.01, 1.0), 0.02;
    replaced.R << 0.25;
    replaced.P0 << 1.5, -0.2, std::nextafter(-0.2, 0.0), 0.8;

    Eigen::Matrix2d drifted_Q = replaced.Q;
    Eigen::Matrix2d drifted_P0 = replaced.P0;
    replaced.Q(0, 1) = replaced.Q(1, 0) = drifted_Q(0, 1) + (drifted_Q(1, 0) - drifted_Q(0, 1)) / 2.0;
    replaced.P0(0, 1) = replaced.P0(1, 0) = drifted_P0(0, 1) + (drifted_P0(1, 0) - drifted_P0(0, 1)) / 2.0;

    auto filter = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, cfg));
    auto twin = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, replaced));

    REQUIRE(filter.set_noise(drifted_Q, replaced.R).has_value());
    REQUIRE(filter.reset_covariance(drifted_P0).has_value());
    CHECK((filter.covariance().array() == replaced.P0.array()).all());

    const Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();
    const Eigen::Matrix<double, 1, 1> z = Eigen::Matrix<double, 1, 1>::Constant(1.0);
    for(int k = 0; k < 5; ++k)
    {
        filter.predict(u);
        twin.predict(u);
        REQUIRE(filter.update(z).has_value());
        REQUIRE(twin.update(z).has_value());
    }
    CHECK((filter.state().array() == twin.state().array()).all());
    CHECK((filter.covariance().array() == twin.covariance().array()).all());
}

TEST_CASE("Kalman rejects correction overflow without committing it",
          "[kalman][hardening][negative]")
{
    auto sys = make_const_velocity_system();
    ctrlpp::kalman_config<double, 2, 1, 1> config;
    config.x0 << -std::numeric_limits<double>::max(), 0.0;
    auto filter = ctrlpp::test::constructed(
        ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, config));
    auto const state_before = filter.state();
    auto const covariance_before = filter.covariance();
    auto const innovation_before = filter.innovation();
    auto const nis_before = filter.nis();

    Eigen::Matrix<double, 1, 1> measurement;
    measurement << std::numeric_limits<double>::max();
    auto result = filter.update(measurement);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::kalman_update_error::non_finite_result);
    CHECK(filter.state() == state_before);
    CHECK(filter.covariance() == covariance_before);
    CHECK(filter.innovation() == innovation_before);
    CHECK(filter.nis() == nis_before);
    CHECK(filter.health() == ctrlpp::kalman_health::ok);
}

TEST_CASE("Kalman zero measurement noise gives the deadbeat gain",
          "[kalman][hardening][negative]")
{
    auto sys = make_const_velocity_system();
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    Eigen::Matrix<double, 1, 1> R;
    R << 0.0;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Identity() * 10.0;

    auto kf = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}));

    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();
    kf.predict(u);

    Eigen::Matrix<double, 1, 1> z;
    z << 5.0;
    REQUIRE(kf.update(z).has_value());

    // With R = 0 the innovation covariance is C P C', still invertible, so the
    // gain in the observed coordinate is exactly one: the measurement is
    // believed completely. Both assertions below are BITWISE, and both are
    // structural rather than lucky.
    //
    // The estimate: the prior estimate in the observed coordinate is exactly
    // zero -- x0 is zero, the input is zero and the state matrix is unit
    // diagonal -- so the correction 0 + 1 * (5 - 0) is the measurement itself
    // with no cancellation anywhere. Finiteness, which this replaces, held for
    // every possible gain including a zero one.
    REQUIRE(kf.state()[0] == 5.0);

    // The covariance: the Joseph factor I - K C has an exactly zero first row,
    // because the gain's observed entry is exactly one, so the first row and
    // column of the posterior are exact zeros and the measurement-noise term
    // K R K' is an exact zero as well.
    REQUIRE(kf.covariance()(0, 0) == 0.0);

    // The symmetrization the filter applies makes the off-diagonal entries the
    // same double, so symmetry is an equality and not a tolerance.
    REQUIRE(kf.covariance()(0, 1) == kf.covariance()(1, 0));
}

TEST_CASE("Kalman NaN measurement is rejected without touching the estimate",
          "[kalman][hardening][negative]")
{
    auto sys = make_const_velocity_system();
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Identity();

    auto kf = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}));
    // Stepped only with the valid measurement, never with the poisoned one, so
    // it says what the filter would have carried had the bad sample never
    // arrived.
    auto reference = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}));

    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();
    kf.predict(u);
    reference.predict(u);

    // Snapshot immediately before the poisoned step.
    const Eigen::Vector2d x_before = kf.state();
    const Eigen::Matrix<double, 2, 2> P_before = kf.covariance();

    Eigen::Matrix<double, 1, 1> z_bad;
    z_bad << std::numeric_limits<double>::quiet_NaN();

    const auto rejected = kf.update(z_bad);

    REQUIRE_FALSE(rejected.has_value());
    REQUIRE(rejected.error() == ctrlpp::kalman_update_error::non_finite_measurement);

    // Exact comparison, not a tolerance: a rejected step performs no arithmetic
    // on the carried estimate at all, so bitwise equality is the contract and a
    // tolerance would admit a step that partially ran.
    CHECK(kf.state() == x_before);
    // The covariance was already measurement-independent before the guard
    // existed -- update_covariance(K) takes only the gain, never z -- so this
    // half of the invariant is structural. The state half is what the guard
    // adds.
    CHECK(kf.covariance() == P_before);
    // A rejection describes the sample, not the filter: nothing was mutated, so
    // the filter is not degraded and must not report that it is.
    CHECK(kf.health() == ctrlpp::kalman_health::ok);

    // The poison did not latch: the next valid step produces exactly what it
    // would have produced had the poisoned step never been attempted.
    Eigen::Matrix<double, 1, 1> z_good;
    z_good << 5.0;
    REQUIRE(kf.update(z_good).has_value());
    REQUIRE(reference.update(z_good).has_value());

    CHECK(kf.state() == reference.state());
    CHECK(kf.covariance() == reference.covariance());
}

TEST_CASE("Kalman zero Q contracts the covariance volume by exactly R over S",
          "[kalman][hardening][negative]")
{
    auto sys = make_const_velocity_system();
    auto Q = Eigen::Matrix<double, 2, 2>::Zero();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Identity();

    auto kf = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}));

    // Two exact statements hold with a zero process noise, and both are asserted
    // at EVERY step rather than at the end:
    //
    //  1. The prediction is a pure congruence, P_pred = A P A', and the state
    //     matrix here has unit determinant, so det(P_pred) = det(P_post_prev).
    //     The optimal update then multiplies the determinant by det(I - K C),
    //     which for a single output is exactly R/S. Together:
    //         det(P_post) * S = det(P_pred) * R.
    //  2. The update never increases the covariance: P_pred - P_post is positive
    //     semidefinite.
    //
    // The estimate is deliberately NOT asserted to reach the constant
    // measurement. It does not: the covariance falls toward zero, the gain with
    // it, and the filter stops listening while the velocity coordinate still
    // carries a bias. After a hundred steps the position estimate is 5.096, not
    // 5, and an oracle claiming convergence to 5 would be false.
    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();
    for(int k = 0; k < 100; ++k)
    {
        CAPTURE(k);
        kf.predict(u);
        const Eigen::Matrix2d P_pred = kf.covariance();
        const auto step = gain_at(sys, P_pred, R(0, 0));

        Eigen::Matrix<double, 1, 1> z;
        z << 5.0;
        REQUIRE(kf.update(z).has_value());
        const Eigen::Matrix2d P_post = kf.covariance();

        // The scale of each determinant is the sum of the magnitudes of the two
        // products it differences, because that difference is where the
        // cancellation lives; a budget written against the determinant itself
        // would be measuring the wrong quantity.
        const double cancel_pred =
            std::abs(P_pred(0, 0) * P_pred(1, 1)) + std::abs(P_pred(0, 1) * P_pred(1, 0));
        const double cancel_post =
            std::abs(P_post(0, 0) * P_post(1, 1)) + std::abs(P_post(0, 1) * P_post(1, 0));
        const double budget =
            kf_determinant_ops * kf_eps
            * std::max(cancel_post * step.innovation_covariance, cancel_pred * R(0, 0));

        REQUIRE_THAT(P_post.determinant() * step.innovation_covariance,
                     WithinAbs(P_pred.determinant() * R(0, 0), budget));

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigsolver(P_pred - P_post);
        REQUIRE(eigsolver.eigenvalues().minCoeff() >= -kf_eig_ops * kf_eps * P_pred.norm());

        REQUIRE(P_post(0, 1) == P_post(1, 0));
    }
}

TEST_CASE("Kalman scalar analytical gain comparison", "[kalman][hardening][precision]")
{
    // Scalar system: A=1, B=0, C=1, D=0 (random walk)
    ctrlpp::discrete_state_space<double, 1, 1, 1> sys;
    sys.A << 1.0;
    sys.B << 0.0;
    sys.C << 1.0;
    sys.D << 0.0;

    double q_val = 0.1;
    double r_val = 1.0;
    double p0 = 10.0;

    Eigen::Matrix<double, 1, 1> Q, R, P0;
    Q << q_val;
    R << r_val;
    P0 << p0;
    Eigen::Matrix<double, 1, 1> x0;
    x0 << 0.0;

    auto kf = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 1, 1, 1>::create(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}));

    // Predict step: P_pred = A*P*A^T + Q = P + Q
    // After first predict: P_pred = p0 + q = 10.1
    // K = P_pred * H' / (H * P_pred * H' + R) = 10.1 / (10.1 + 1.0) = 10.1/11.1
    Eigen::Matrix<double, 1, 1> u_zero;
    u_zero << 0.0;
    kf.predict(u_zero);

    double P_pred = p0 + q_val;
    double expected_K = P_pred / (P_pred + r_val);

    Eigen::Matrix<double, 1, 1> z;
    z << 5.0;
    REQUIRE(kf.update(z).has_value());

    // After update: x = 0 + K*(5 - 0) = K*5.
    //
    // The derivation above is the strongest oracle in this file and only its
    // tolerance was unprincipled. Four rounded operations reach the reference
    // value -- the sum forming the predicted covariance, the sum forming the
    // innovation denominator, the division, and the product against the
    // measurement -- and the filter reaches the same value through a one-by-one
    // QR solve whose backward error is bounded by two more.
    constexpr int scalar_update_ops = 6;
    REQUIRE_THAT(kf.state()[0], WithinRel(expected_K * 5.0, scalar_update_ops * kf_eps));
}

TEST_CASE("Kalman covariance stays positive definite over 1000 steps",
          "[kalman][hardening][stability]")
{
    auto sys = make_const_velocity_system();
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Identity() * 10.0;

    auto kf = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}));

    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();

    for(int k = 0; k < 1000; ++k)
    {
        CAPTURE(k);
        kf.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << 5.0 + 0.01 * k;
        REQUIRE(kf.update(z).has_value());

        // Symmetry is the other half of the covariance contract and was checked
        // nowhere. It is an EQUALITY here rather than a tolerance: the filter
        // symmetrizes explicitly, so the two off-diagonal entries are formed as
        // the same expression and are the same double.
        REQUIRE(kf.covariance()(0, 1) == kf.covariance()(1, 0));

        // The floor is the symmetric eigensolver's own backward error at the
        // scale of the matrix it was handed, not a round number below zero. The
        // threshold this replaces, -1e-10, is about 450000 epsilons of slack in
        // the WRONG direction: it admitted a genuinely indefinite covariance for
        // a thousand consecutive steps, which is the opposite of what the case
        // is named for.
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> eigsolver(kf.covariance());
        REQUIRE(eigsolver.eigenvalues().minCoeff() >= -kf_eig_ops * kf_eps * kf.covariance().norm());
    }
}

TEST_CASE("Kalman estimation error follows its own closed-loop recursion exactly",
          "[kalman][hardening][convergence]")
{
    auto sys = make_const_velocity_system();
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    Eigen::Matrix<double, 1, 1> R;
    R << 0.1;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Identity() * 10.0;

    auto kf = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}));

    // The truth is a full state, advanced by the same matrix the filter models,
    // and every measurement is that state's observed coordinate with no noise.
    // The estimation error is therefore homogeneous:
    //     e(k+1) = (I - K(k) C) A e(k),   e(0) = x_hat(0) - x_true(0),
    // and the test propagates it alongside using the gain implied by the
    // filter's own reported predicted covariance.
    Eigen::Vector2d x_true;
    x_true << 0.0, 1.0;
    Eigen::Vector2d error = x0 - x_true;

    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();

    // The realized error and the propagated one are two floating-point
    // realizations of the same recursion, so their divergence is injected once
    // per step at the counted step rounding and then carried forward by the same
    // closed-loop map. The bound is accumulated with that map's own spectral
    // norm rather than assumed non-amplifying: it exceeds one during the
    // transient, so a plain multiple of the step count would not bound it.
    double divergence_budget = 0.0;

    for(int k = 0; k < 200; ++k)
    {
        CAPTURE(k);
        x_true = (sys.A * x_true).eval();
        kf.predict(u);

        const auto step = gain_at(sys, kf.covariance(), R(0, 0));
        error = (step.closed_loop * error).eval();

        Eigen::Matrix<double, 1, 1> z;
        z << x_true(0);
        REQUIRE(kf.update(z).has_value());

        const double scale = std::max({x_true.cwiseAbs().maxCoeff(),
                                       kf.state().cwiseAbs().maxCoeff(), 1.0});
        divergence_budget = step.closed_loop.norm() * divergence_budget + kf_step_ops * kf_eps * scale;

        const Eigen::Vector2d realized = kf.state() - x_true;
        CAPTURE(realized(0), realized(1), error(0), error(1));
        // Two-sided by construction, which the threshold this replaces -- a bare
        // 0.5 against a realized error of 7e-12 -- could not be: a filter that
        // converged faster than its own gain sequence permits fails this too.
        REQUIRE_THAT(realized(0), WithinAbs(error(0), divergence_budget));
        REQUIRE_THAT(realized(1), WithinAbs(error(1), divergence_budget));
    }
}

TEST_CASE("Kalman unobservable mode diverges linearly while the observed one converges",
          "[kalman][hardening][robustness]")
{
    // The framing this case used to carry -- ill-conditioning -- named the wrong
    // property. The builder gives A = diag(a, 1) and C = [1, 0], so the SECOND
    // mode has eigenvalue exactly one and is not observed at all. Its covariance
    // is driven by the process noise and by nothing else, so it grows by exactly
    // Q(1,1) every step and diverges linearly. Finiteness hid a covariance that
    // is provably diverging, and would have kept hiding it for any number of
    // steps.
    constexpr double a = 1e-10;
    constexpr int steps = 100;
    auto sys = ctrlpp::test::near_singular_system<double, 2, 1, 1>(a);
    Eigen::Matrix<double, 2, 2> Q = Eigen::Matrix<double, 2, 2>::Identity() * 0.01;
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 2, 2> P0 = Eigen::Matrix<double, 2, 2>::Identity();

    auto kf = ctrlpp::test::constructed(ctrlpp::kalman_filter<double, 2, 1, 1>::create(sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}));

    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();

    for(int k = 0; k < steps; ++k)
    {
        kf.predict(u);
        Eigen::Matrix<double, 1, 1> z;
        z << 1.0;
        REQUIRE(kf.update(z).has_value());
    }

    // The unobservable coordinate is never corrected, and it starts at exactly
    // zero: bitwise.
    REQUIRE(kf.state()[1] == 0.0);

    // The two modes never couple, because A, C, Q and P0 are all diagonal and
    // every product in the recursion preserves that: bitwise zero off-diagonals.
    REQUIRE(kf.covariance()(0, 1) == 0.0);
    REQUIRE(kf.covariance()(1, 0) == 0.0);

    // Linear growth in closed form: one addition of the process-noise entry per
    // step, on an accumulator whose largest value is the answer itself.
    constexpr int growth_ops = steps;
    REQUIRE_THAT(kf.covariance()(1, 1),
                 WithinRel(P0(1, 1) + steps * Q(1, 1), growth_ops * kf_eps));

    // The observed coordinate reaches a closed-form fixed point. The predicted
    // covariance there settles at Q(0,0), because the congruence by a shrinks
    // the carried posterior below one ulp of it, so the gain settles at
    // K = Q(0,0) / (Q(0,0) + R); the estimate then satisfies
    // x = a*x*(1 - K) + K, whose solution is written below. The steps needed to
    // reach it are two, since a squared is already negligible.
    const double P_pred_00 = Q(0, 0);
    const double K = P_pred_00 / (P_pred_00 + R(0, 0));
    const double fixed_point = K / (1.0 - a * (1.0 - K));
    //
    // Note what that fixed point is NOT: a filter told 1.0 a hundred times
    // reports 0.0099, because the state matrix annihilates the estimate between
    // steps and only the process noise keeps the gain from vanishing. Asserting
    // the closed form pins that; asserting convergence to the measurement would
    // have been false.
    constexpr int fixed_point_ops = 8;
    REQUIRE_THAT(kf.state()[0], WithinRel(fixed_point, fixed_point_ops * kf_eps));
}

TEST_CASE("Kalman is_steady_state with near-zero covariance",
          "[kalman][hardening][coverage]")
{
    constexpr std::size_t NX = 2, NU = 1, NY = 1;
    using KF = ctrlpp::kalman_filter<double, NX, NU, NY>;

    Eigen::Matrix2d A;
    A << 1.0, 0.01, 0.0, 1.0;
    Eigen::Vector2d B(0.0, 0.01);
    Eigen::RowVector2d C(1.0, 0.0);
    Eigen::Matrix<double, 1, 1> D = Eigen::Matrix<double, 1, 1>::Zero();
    ctrlpp::discrete_state_space<double, NX, NU, NY> sys{A, B, C, D};

    ctrlpp::kalman_config<double, NX, NU, NY> cfg{};
    // Near-zero initial covariance triggers P.norm() < epsilon branch
    cfg.P0 = Eigen::Matrix2d::Zero();
    cfg.Q = Eigen::Matrix2d::Identity() * 1e-300;
    cfg.R = Eigen::Matrix<double, 1, 1>::Identity();

    auto kf = ctrlpp::test::constructed(KF::create(sys, cfg));

    // With zero P, should immediately report steady state
    CHECK(kf.is_steady_state());

    // Pinned from the other side as well, because a predicate that answers yes
    // to everything is not a predicate. A large initial covariance driven by a
    // large process noise is emphatically NOT at steady state after one step:
    // the posterior differs from the prior posterior by a relative amount many
    // orders above the predicate's own tolerance.
    ctrlpp::kalman_config<double, NX, NU, NY> moving{};
    moving.P0 = Eigen::Matrix2d::Identity() * 1e6;
    moving.Q = Eigen::Matrix2d::Identity();
    moving.R = Eigen::Matrix<double, 1, 1>::Identity();

    auto transient = ctrlpp::test::constructed(KF::create(sys, moving));
    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();
    transient.predict(u);
    Eigen::Matrix<double, 1, 1> z;
    z << 1.0;
    REQUIRE(transient.update(z).has_value());
    CHECK_FALSE(transient.is_steady_state());
}
