#include "hardening_helpers.h"
#include "ctrlpp/control/dare.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;

TEST_CASE("DARE NaN in A returns nullopt or NaN", "[dare][hardening][negative]")
{
    auto A = ctrlpp::test::nan_matrix<double, 2, 2>();
    Eigen::Matrix<double, 2, 1> B;
    B << 1.0, 0.0;
    auto Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    CHECK(!result.has_value());
}

TEST_CASE("DARE NaN in B returns nullopt", "[dare][hardening][negative]")
{
    auto A = Eigen::Matrix<double, 2, 2>::Identity();
    auto B = ctrlpp::test::nan_matrix<double, 2, 1>();
    auto Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    CHECK(!result.has_value());
}

TEST_CASE("DARE refuses a singular R", "[dare][hardening][negative]")
{
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 1.0, 0.0, 1.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.5, 1.0;
    auto Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 0.0;

    auto const result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE_FALSE(result.has_value());

    // The enumerator is non_finite_input even though every argument is finite,
    // and that is what the enumerator's own definition says: it covers "A, B, Q,
    // R OR the assembled symplectic Z contains NaN/Inf". The symplectic build
    // needs R^{-1} to form G = B R^{-1} B', and a singular R makes that whole
    // block infinite -- Eigen's rank-revealing QR solve of the 1x1 zero returns
    // infinity rather than declining -- so Z fails its finiteness test.
    //
    // The refusal is right. Note what it does NOT say: the module has
    // singular_a for the analogous condition on the state matrix, and no
    // counterpart naming a singular weighting matrix, so the caller is told to
    // look for a non-finite input when the actual obstacle is a weighting they
    // set to zero on purpose.
    CHECK(result.error() == ctrlpp::dare_error::non_finite_input);
}

TEST_CASE("DARE known 2x2 solution is positive definite", "[dare][hardening][precision]")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    A << 1.0, 1.0, 0.0, 1.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.5, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    auto& P = result->P;

    // Verify positive definite
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> eigsolver(P);
    for(int i = 0; i < 2; ++i)
        CHECK(eigsolver.eigenvalues()(i) > 0.0);

    // Verify symmetric
    REQUIRE_THAT((P - P.transpose()).norm(), WithinAbs(0.0, 1e-14));
}

TEST_CASE("DARE scalar analytical solution", "[dare][hardening][precision]")
{
    Eigen::Matrix<double, 1, 1> A, B, Q, R;
    A(0, 0) = 1.0;
    B(0, 0) = 1.0;
    Q(0, 0) = 1.0;
    R(0, 0) = 1.0;

    auto result = ctrlpp::dare<double, 1, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    // Analytical: P = golden ratio = (1+sqrt(5))/2
    double golden = (1.0 + std::sqrt(5.0)) / 2.0;
    REQUIRE_THAT(result->P(0, 0), WithinAbs(golden, 1e-10));
}

TEST_CASE("DARE solution is positive definite for stable system", "[dare][hardening][stability]")
{
    Eigen::Matrix<double, 2, 2> A, Q;
    A << 0.9, 0.1, 0.0, 0.8;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.0, 1.0;
    Q = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> eigsolver(result->P);
    for(int i = 0; i < 2; ++i)
        CHECK(eigsolver.eigenvalues()(i) > -1e-10);
}

TEST_CASE("DARE solves an ill-conditioned but well-posed problem",
          "[dare][hardening][robustness]")
{
    // A = [[1,1],[0,1]] is controllable from B = [0.5; 1] (rank[B, AB] = 2), and
    // Q = diag(1, 1e-10) is positive definite -- barely -- which makes the pair
    // detectable. A unique stabilizing positive-definite solution therefore
    // exists, so a refusal here would be a well-posed problem reported
    // unsolvable, and finiteness alone would be far weaker than the property
    // the solution is supposed to have.
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 1.0, 0.0, 1.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.5, 1.0;
    auto Q = ctrlpp::test::ill_conditioned_2x2<double>(1e10);
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto const result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    REQUIRE(result.has_value());

    auto const& P = result->P;
    constexpr double eps = std::numeric_limits<double>::epsilon();

    auto const res = ctrlpp::test::riccati_residual<double, 2, 1>(A, B, Q, R, P);
    CAPTURE(res.norm, res.scale);
    REQUIRE(res.norm <= ctrlpp::test::riccati_residual_ops<2, 1> * eps * res.scale);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> pes(P);
    for(int i = 0; i < 2; ++i)
        CHECK(pes.eigenvalues()(i) > 0.0);

    auto const K = ctrlpp::test::riccati_gain<double, 2, 1>(A, B, R, P);
    Eigen::Matrix<double, 2, 2> Acl = (A - B * K).eval();
    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2>> ces(Acl, false);
    for(int i = 0; i < 2; ++i)
        REQUIRE(std::abs(ces.eigenvalues()(i)) < 1.0);
}

TEST_CASE("DARE NaN in Q returns nullopt", "[dare][hardening][negative]")
{
    Eigen::Matrix<double, 2, 2> A;
    A << 1.0, 1.0, 0.0, 1.0;
    Eigen::Matrix<double, 2, 1> B;
    B << 0.5, 1.0;
    auto Q = ctrlpp::test::nan_matrix<double, 2, 2>();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto result = ctrlpp::dare<double, 2, 1>(A, B, Q, R);
    CHECK(!result.has_value());
}
