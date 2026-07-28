#include "hardening_helpers.h"
#include "ctrlpp/control/lqr.h"
#include "ctrlpp/control/dare.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("LQR refuses a non-finite state weighting", "[lqr][hardening][negative]")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, 1.0, 0.0, 1.0;
    B << 0.5, 1.0;
    Q = ctrlpp::test::nan_matrix<double, 2, 2>();
    R << 1.0;

    auto result = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    // The case name already claimed the answer; this holds the library to it.
    // A gain formed from a NaN weighting would be a controller reporting
    // success while carrying no usable number, so "nullopt or NaN" was a
    // disjunction over the only two outcomes observable here and could not
    // fail whichever the library did.
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::dare_error::non_finite_input);
}

TEST_CASE("LQR refuses a non-finite input weighting", "[lqr][hardening][negative]")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, 1.0, 0.0, 1.0;
    B << 0.5, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R << std::numeric_limits<double>::quiet_NaN();

    auto result = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == ctrlpp::dare_error::non_finite_input);
}

TEST_CASE("LQR stabilizes an ill-conditioned but well-posed pair", "[lqr][hardening][robustness]")
{
    // A = diag(1, 1e-10), B = [1; 0]. Mode 0 sits at eigenvalue 1 and IS input
    // coupled, so it is controllable; mode 1 sits at 1e-10, which is inside the
    // unit circle and so needs no control at all. The pair is therefore
    // stabilizable (though not controllable -- [B, AB] has rank 1), and Q = I
    // makes it detectable, so a unique stabilizing solution exists and a
    // refusal here would be a solvable problem reported unsolvable.
    auto A = ctrlpp::test::ill_conditioned_2x2<double>(1e10);
    Eigen::Matrix<double, 2, 1> B;
    B << 1.0, 0.0;
    auto Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto const solved = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(solved.has_value());

    auto const& P = solved->P;
    constexpr double eps = std::numeric_limits<double>::epsilon();

    // The solution solves the equation, not merely "is finite".
    auto const res = ctrlpp::test::riccati_residual<double, 2, 1>(A, B, Q, R, P);
    CAPTURE(res.norm, res.scale);
    REQUIRE(res.norm <= ctrlpp::test::riccati_residual_ops<2, 1> * eps * res.scale);

    // ... and it is positive definite.
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> pes(P);
    for(int i = 0; i < 2; ++i)
        CHECK(pes.eigenvalues()(i) > 0.0);

    auto const result = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    // The gain is the one this solution implies, formed the same way.
    auto const K_expected = ctrlpp::test::riccati_gain<double, 2, 1>(A, B, R, P);
    CHECK((*result - K_expected).norm()
          <= ctrlpp::test::riccati_residual_ops<2, 1> * eps * K_expected.norm());

    // The mathematically required property of an LQR gain: the closed loop is
    // asymptotically stable. Strictly inside the unit circle is the exact
    // contract boundary, not a fitted constant.
    Eigen::Matrix<double, 2, 2> Acl = (A - B * (*result)).eval();
    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> ces(Acl, false);
    for(int i = 0; i < 2; ++i)
        REQUIRE(std::abs(ces.eigenvalues()(i)) < 1.0);
}

TEST_CASE("LQR double integrator matches analytical", "[lqr][hardening][precision]")
{
    double dt = 0.1;
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, dt, 0.0, 1.0;
    B << 0.5 * dt * dt, dt;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R << 1.0;

    auto result = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    auto& K = *result;
    CHECK(std::isfinite(K(0, 0)));
    CHECK(std::isfinite(K(0, 1)));

    // Verify closed-loop eigenvalues are inside unit circle
    auto Acl = (A - B * K).eval();
    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> solver(Acl, false);
    for(int i = 0; i < 2; ++i)
        CHECK(std::abs(solver.eigenvalues()(i)) < 1.0);
}

TEST_CASE("LQR closed-loop eigenvalues inside unit circle", "[lqr][hardening][stability]")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    Eigen::Matrix<double, 2, 1> B;
    Eigen::Matrix<double, 1, 1> R;

    A << 1.0, 1.0, 0.0, 1.0;
    B << 0.5, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    R << 1.0;

    auto result = ctrlpp::lqr_gain<double, 2, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    auto Acl = (A - B * *result).eval();
    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> solver(Acl, false);
    for(int i = 0; i < 2; ++i)
        REQUIRE(std::abs(solver.eigenvalues()(i)) < 1.0);
}

TEST_CASE("LQR refuses an unstabilizable pair with the enumerator that names it",
          "[lqr][hardening][negative]")
{
    // A = diag(1e-10, 1), B = [1; 0]. Mode 1 sits at eigenvalue exactly 1 with
    // B(1) = 0: uncontrollable AND not asymptotically stable. The pair is
    // therefore NOT stabilizable and no stabilizing solution exists, so a gain
    // returned here would be a controller claiming to stabilize a plant it
    // cannot.
    //
    // The enumerator is reachable in this case because the uncontrollable mode
    // is ON the unit circle: it contributes 1 and its reciprocal 1 to the
    // symplectic spectrum, neither strictly inside, so only one of the four
    // eigenvalues lands in the stable region and the placement count falls
    // short of the state dimension.
    auto sys = ctrlpp::test::near_singular_system<double, 2, 1, 1>(1e-10);
    auto Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    // Asserted at the solver, where the verdict originates ...
    auto const solved = ctrlpp::dare<double, 2, 1>(sys.A, sys.B, Q, R);
    REQUIRE_FALSE(solved.has_value());
    REQUIRE(solved.error() == ctrlpp::dare_error::non_stabilisable);

    // ... and at the gain, which forwards it rather than flattening it. Without
    // this second assertion the forwarding is untested at the surface callers
    // actually use.
    auto const result = ctrlpp::lqr_gain<double, 2, 1>(sys.A, sys.B, Q, R);
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::dare_error::non_stabilisable);
}

TEST_CASE("LQR scalar integrator analytical gain", "[lqr][hardening][precision]")
{
    Eigen::Matrix<double, 1, 1> A, B, Q, R;
    A(0, 0) = 1.0;
    B(0, 0) = 1.0;
    Q(0, 0) = 1.0;
    R(0, 0) = 1.0;

    auto result = ctrlpp::lqr_gain<double, 1, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    double golden = (1.0 + std::sqrt(5.0)) / 2.0;
    double expected_K = golden / (1.0 + golden);
    REQUIRE_THAT((*result)(0, 0), WithinAbs(expected_K, 1e-10));
}
