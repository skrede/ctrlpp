#include "ctrlpp/detail/riccati_solution.h"


#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/QR>
#include <Eigen/Dense>

#include <cmath>
#include <limits>

using Catch::Matchers::WithinAbs;


TEST_CASE("extract_riccati_solution recovers symmetric PSD P on well-conditioned input",
          "[riccati_solution][success]")
{
    constexpr int n  = 2;
    constexpr int N2 = 2 * n;
    using Scalar = double;

    // Build a real-orthogonal basis whose leading n columns span the same
    // invariant subspace as [I; R_target]. After Householder-QR, the full
    // 2n x 2n Q_full satisfies (Q_full.block(n, 0) * Q_full.block(0, 0)^-1)
    // = R_target (same column span, same ratio).
    Eigen::Matrix<Scalar, n, n> R_target;
    R_target << 2.0, 0.5, 0.5, 3.0;

    Eigen::Matrix<Scalar, N2, n> M;
    M.template block<n, n>(0, 0) = Eigen::Matrix<Scalar, n, n>::Identity();
    M.template block<n, n>(n, 0) = R_target;

    Eigen::HouseholderQR<Eigen::Matrix<Scalar, N2, n>> qr(M);
    Eigen::Matrix<Scalar, N2, N2> Q_full = qr.householderQ();

    auto result = ctrlpp::detail::extract_riccati_solution<Scalar, N2>(Q_full);
    REQUIRE(result.has_value());

    const auto& P = *result;
    const Scalar tol =
        100.0 * std::numeric_limits<Scalar>::epsilon() * R_target.cwiseAbs().maxCoeff();
    CHECK((P - R_target).norm() < tol);
    CHECK((P - P.transpose()).norm() < tol);
}

TEST_CASE("extract_riccati_solution returns singular_u11 when U11 is rank-deficient",
          "[riccati_solution][error]")
{
    constexpr int N2 = 4;
    using Scalar = double;

    // U11 is the leading 2x2 block; make it identically zero while keeping
    // the full 4x4 U orthogonal by permuting basis vectors.
    Eigen::Matrix<Scalar, N2, N2> U = Eigen::Matrix<Scalar, N2, N2>::Zero();
    U(2, 0) = 1.0;
    U(3, 1) = 1.0;
    U(0, 2) = 1.0;
    U(1, 3) = 1.0;

    auto result = ctrlpp::detail::extract_riccati_solution<Scalar, N2>(U);
    REQUIRE(!result.has_value());
    CHECK(result.error() == ctrlpp::detail::riccati_extract_error::singular_u11);
}

TEST_CASE("extract_riccati_solution returns non_psd for negative-definite P",
          "[riccati_solution][error]")
{
    constexpr int N2 = 4;
    using Scalar = double;

    // U11 = I, U21 = -I  ->  P = -I, which is negative definite.
    Eigen::Matrix<Scalar, N2, N2> U = Eigen::Matrix<Scalar, N2, N2>::Zero();
    U(0, 0) = 1.0;
    U(1, 1) = 1.0;
    U(2, 0) = -1.0;
    U(3, 1) = -1.0;
    U(2, 2) = 1.0;
    U(3, 3) = 1.0;

    auto result = ctrlpp::detail::extract_riccati_solution<Scalar, N2>(U);
    REQUIRE(!result.has_value());
    CHECK(result.error() == ctrlpp::detail::riccati_extract_error::non_psd);
}

TEST_CASE("extract_riccati_solution handles NaN in U via non_finite or non_psd",
          "[riccati_solution][error]")
{
    constexpr int N2 = 4;
    using Scalar = double;

    Eigen::Matrix<Scalar, N2, N2> U = Eigen::Matrix<Scalar, N2, N2>::Identity();
    U(2, 0) = std::numeric_limits<Scalar>::quiet_NaN();

    auto result = ctrlpp::detail::extract_riccati_solution<Scalar, N2>(U);
    REQUIRE(!result.has_value());
    CHECK((result.error() == ctrlpp::detail::riccati_extract_error::non_finite
        || result.error() == ctrlpp::detail::riccati_extract_error::non_psd
        || result.error() == ctrlpp::detail::riccati_extract_error::singular_u11));
}
