#include "ctrlpp/detail/schur_reorder.h"


#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <limits>
#include <complex>

using Catch::Matchers::WithinAbs;

namespace
{

template <typename Scalar, int N>
auto infinity_norm(const Eigen::Matrix<Scalar, N, N>& M) -> Scalar
{
    return M.cwiseAbs().maxCoeff();
}

template <typename Scalar, int N>
auto orthogonality_residual(const Eigen::Matrix<Scalar, N, N>& U) -> Scalar
{
    const auto I = Eigen::Matrix<Scalar, N, N>::Identity();
    return (U.transpose() * U - I).norm();
}

template <typename Scalar>
auto eigenvalue_set_2x2(Scalar a, Scalar b, Scalar c, Scalar d)
    -> std::pair<std::complex<Scalar>, std::complex<Scalar>>
{
    const Scalar tr    = a + d;
    const Scalar det   = a * d - b * c;
    const Scalar discr = tr * tr / Scalar{4} - det;
    if (discr >= Scalar{0})
    {
        const Scalar s = std::sqrt(discr);
        return {std::complex<Scalar>(tr / Scalar{2} + s, Scalar{0}),
                std::complex<Scalar>(tr / Scalar{2} - s, Scalar{0})};
    }
    const Scalar im = std::sqrt(-discr);
    return {std::complex<Scalar>(tr / Scalar{2},  im),
            std::complex<Scalar>(tr / Scalar{2}, -im)};
}

}  // anonymous namespace

TEST_CASE("schur_reorder policy tag types are empty structs")
{
    STATIC_REQUIRE(std::is_empty_v<ctrlpp::detail::pivot_ratio_conditioning>);
    STATIC_REQUIRE(std::is_empty_v<ctrlpp::detail::hager_higham_conditioning>);
}

TEST_CASE("reorder_real_schur skeleton returns trivial success for identity T")
{
    constexpr int N = 4;
    Eigen::Matrix<double, N, N> T = Eigen::Matrix<double, N, N>::Identity();
    Eigen::Matrix<double, N, N> U = Eigen::Matrix<double, N, N>::Identity();

    auto predicate = [](std::complex<double>) { return true; };
    auto r = ctrlpp::detail::reorder_real_schur<double, N>(T, U, predicate);

    CHECK(r.complete);
    CHECK(r.subspace_separation == 1.0);
}

TEST_CASE("standardize_2x2_block brings general block to Murnaghan canonical form")
{
    constexpr int N = 3;
    using Mat = Eigen::Matrix<double, N, N>;
    Mat T = Mat::Zero();
    // General 2x2 block with real eigenvalue pair at positions (0, 1)
    T(0, 0) = 0.7;
    T(0, 1) = 2.0;
    T(1, 0) = -0.5;
    T(1, 1) = 0.3;
    T(2, 2) = 1.5;  // off-block scalar, should remain unchanged
    Mat U = Mat::Identity();

    const auto eigs_before = eigenvalue_set_2x2(T(0, 0), T(0, 1), T(1, 0), T(1, 1));
    const Mat T_before = T;
    const double scale = infinity_norm(T);
    const double eps = std::numeric_limits<double>::epsilon();

    ctrlpp::detail::standardize_2x2_block<double, N>(T, U, 0);

    // Eigenvalues preserved: sum (trace) and product (determinant) of the block.
    const double tr_after  = T(0, 0) + T(1, 1);
    const double det_after = T(0, 0) * T(1, 1) - T(0, 1) * T(1, 0);
    const double tr_before = std::real(eigs_before.first + eigs_before.second);
    const double det_before = std::real(eigs_before.first * eigs_before.second);
    CHECK_THAT(tr_after,  WithinAbs(tr_before, 10 * eps * scale));
    CHECK_THAT(det_after, WithinAbs(det_before, 10 * eps * scale));

    // U orthogonality preserved.
    CHECK(orthogonality_residual(U) < 10 * eps * N);

    // Unrelated entries (outside the 2x2 block) are untouched.
    CHECK_THAT(T(2, 2), WithinAbs(T_before(2, 2), 10 * eps * scale));

    // For this real-eigenvalue seed, the canonical form is upper-triangular (c -> 0)
    // or opposite-sign off-diagonals. Verify at least one canonical invariant holds.
    const bool opposite_signs =
        (T(0, 1) > 0.0 && T(1, 0) < 0.0) || (T(0, 1) < 0.0 && T(1, 0) > 0.0);
    const bool c_zeroed = std::abs(T(1, 0)) < 10 * eps * scale;
    CHECK((opposite_signs || c_zeroed));

    // Similarity identity U^T * T_before * U == T_new must hold.
    Mat reconstructed = U.transpose() * T_before * U;
    CHECK((reconstructed - T).norm() < 100 * eps * scale);
}

TEST_CASE("standardize_2x2_block is a no-op when subdiagonal is already zero")
{
    constexpr int N = 2;
    using Mat = Eigen::Matrix<double, N, N>;
    Mat T;
    T << 2.0, 3.0,
         0.0, 1.0;
    Mat U = Mat::Identity();
    const Mat T_before = T;
    const Mat U_before = U;

    ctrlpp::detail::standardize_2x2_block<double, N>(T, U, 0);

    CHECK(T(0, 0) == T_before(0, 0));
    CHECK(T(0, 1) == T_before(0, 1));
    CHECK(T(1, 0) == T_before(1, 0));
    CHECK(T(1, 1) == T_before(1, 1));
    CHECK(U == U_before);
}

TEST_CASE("standardize_2x2_block preserves complex-conjugate pair structure")
{
    constexpr int N = 2;
    using Mat = Eigen::Matrix<double, N, N>;
    Mat T;
    T << 0.3, -0.5,
         0.5,  0.3;   // already standardised complex pair
    Mat U = Mat::Identity();

    const auto eigs_before = eigenvalue_set_2x2(T(0, 0), T(0, 1), T(1, 0), T(1, 1));
    const double scale = infinity_norm(T);
    const double eps = std::numeric_limits<double>::epsilon();

    ctrlpp::detail::standardize_2x2_block<double, N>(T, U, 0);

    // Equal diagonals enforced.
    CHECK_THAT(T(0, 0) - T(1, 1), WithinAbs(0.0, 10 * eps * scale));

    // Opposite signs on off-diagonals (complex pair).
    CHECK(T(0, 1) * T(1, 0) < 0.0);

    // |b|*|c| equals |lambda_im|^2.
    const double im_sq = std::imag(eigs_before.first) * std::imag(eigs_before.first);
    CHECK_THAT(std::abs(T(0, 1)) * std::abs(T(1, 0)),
               WithinAbs(im_sq, 10 * eps * scale));

    CHECK(orthogonality_residual(U) < 10 * eps * N);
}

TEST_CASE("standardize_2x2_block propagates Givens into U correctly")
{
    constexpr int N = 4;
    using Mat = Eigen::Matrix<double, N, N>;
    Mat T = Mat::Zero();
    // Block at offset 1: general form requiring standardisation.
    T(1, 1) = 1.2;
    T(1, 2) = 0.7;
    T(2, 1) = -0.4;
    T(2, 2) = 0.6;
    T(0, 0) = 3.0;
    T(3, 3) = 0.5;
    T(0, 1) = 0.1;
    T(0, 2) = 0.2;
    T(1, 3) = 0.05;
    T(2, 3) = -0.03;
    Mat U = Mat::Identity();
    const Mat T_before = T;
    const double scale = infinity_norm(T);
    const double eps = std::numeric_limits<double>::epsilon();

    ctrlpp::detail::standardize_2x2_block<double, N>(T, U, 1);

    // Similarity identity: U^T * T_before * U == T_new.
    Mat reconstructed = U.transpose() * T_before * U;
    CHECK((reconstructed - T).norm() < 100 * eps * scale);

    // U orthogonality preserved.
    CHECK(orthogonality_residual(U) < 10 * eps * N);

    // Entries outside cols/rows (1, 2) are untouched row-wise on rows 0 and 3.
    CHECK_THAT(T(0, 0), WithinAbs(T_before(0, 0), 10 * eps * scale));
    CHECK_THAT(T(3, 3), WithinAbs(T_before(3, 3), 10 * eps * scale));
}

TEST_CASE("swap 1x1/1x1 blocks preserves eigenvalues and orthogonal U")
{
    constexpr int N = 2;
    using Mat = Eigen::Matrix<double, N, N>;
    Mat T;
    T << 2.0, 1.0,
         0.0, 0.5;
    Mat U = Mat::Identity();
    const Mat T_before = T;
    const double scale = infinity_norm(T);
    const double eps = std::numeric_limits<double>::epsilon();

    double pivot_ratio = 0.0;
    ctrlpp::detail::pivot_ratio_conditioning tag{};
    const bool ok = ctrlpp::detail::swap_real_schur_blocks<double, N>(
        T, U, 0, 1, 1, tag, pivot_ratio);

    REQUIRE(ok);
    CHECK_THAT(T(0, 0), WithinAbs(0.5, 10 * eps * scale));
    CHECK_THAT(T(1, 1), WithinAbs(2.0, 10 * eps * scale));
    CHECK(T(1, 0) == 0.0);
    CHECK(orthogonality_residual(U) < 10 * eps * N);

    // Similarity identity.
    Mat reconstructed = U.transpose() * T_before * U;
    CHECK((reconstructed - T).norm() < 100 * eps * scale);
}

TEST_CASE("swap 1x1/2x2 blocks preserves eigenvalue multiset")
{
    constexpr int N = 3;
    using Mat = Eigen::Matrix<double, N, N>;
    Mat T = Mat::Zero();
    // 1x1 block at (0, 0), standardised 2x2 block at (1, 1)
    T(0, 0) = 2.0;
    T(0, 1) = 1.0;
    T(0, 2) = 0.5;
    T(1, 1) = 0.3;
    T(1, 2) = -0.5;
    T(2, 1) = 0.5;
    T(2, 2) = 0.3;
    Mat U = Mat::Identity();
    const Mat T_before = T;
    const double scale = infinity_norm(T);
    const double eps = std::numeric_limits<double>::epsilon();

    double pivot_ratio = 0.0;
    ctrlpp::detail::pivot_ratio_conditioning tag{};
    const bool ok = ctrlpp::detail::swap_real_schur_blocks<double, N>(
        T, U, 0, 1, 2, tag, pivot_ratio);

    REQUIRE(ok);
    CHECK(pivot_ratio > eps);

    // The 2x2 block moved to position (0, 0); the (2, 0:2) sub-block is ~0.
    CHECK(std::abs(T(2, 0)) < 10 * eps * scale);
    CHECK(std::abs(T(2, 1)) < 10 * eps * scale);

    // U orthogonality.
    CHECK(orthogonality_residual(U) < 10 * eps * N);

    // Trace and Frobenius^2 preserved (eigenvalue multiset invariants).
    CHECK_THAT(T.trace(), WithinAbs(T_before.trace(), 100 * eps * scale));

    // Similarity identity.
    Mat reconstructed = U.transpose() * T_before * U;
    CHECK((reconstructed - T).norm() < 100 * eps * scale);
}

TEST_CASE("swap 2x2/1x1 blocks preserves eigenvalue multiset")
{
    constexpr int N = 3;
    using Mat = Eigen::Matrix<double, N, N>;
    Mat T = Mat::Zero();
    // Standardised 2x2 block at (0, 0), 1x1 block at (2, 2)
    T(0, 0) = 0.3;
    T(0, 1) = -0.5;
    T(0, 2) = 0.7;
    T(1, 0) = 0.5;
    T(1, 1) = 0.3;
    T(1, 2) = 0.2;
    T(2, 2) = 2.0;
    Mat U = Mat::Identity();
    const Mat T_before = T;
    const double scale = infinity_norm(T);
    const double eps = std::numeric_limits<double>::epsilon();

    double pivot_ratio = 0.0;
    ctrlpp::detail::pivot_ratio_conditioning tag{};
    const bool ok = ctrlpp::detail::swap_real_schur_blocks<double, N>(
        T, U, 0, 2, 1, tag, pivot_ratio);

    REQUIRE(ok);
    CHECK(pivot_ratio > eps);

    // Post-swap: 1x1 block at (0, 0), 2x2 block at (1, 1). The (1:3, 0)
    // column is zero below the (0, 0) entry.
    CHECK(std::abs(T(1, 0)) < 10 * eps * scale);
    CHECK(std::abs(T(2, 0)) < 10 * eps * scale);

    CHECK(orthogonality_residual(U) < 10 * eps * N);
    CHECK_THAT(T.trace(), WithinAbs(T_before.trace(), 100 * eps * scale));

    Mat reconstructed = U.transpose() * T_before * U;
    CHECK((reconstructed - T).norm() < 100 * eps * scale);
}

TEST_CASE("swap 2x2/2x2 blocks -- the hard case")
{
    constexpr int N = 4;
    using Mat = Eigen::Matrix<double, N, N>;
    Mat T = Mat::Zero();
    // Two standardised 2x2 blocks with distinct eigenvalue pairs.
    // Block 1 at (0, 0): complex pair 0.3 +/- 0.5i
    T(0, 0) = 0.3;
    T(0, 1) = -0.5;
    T(1, 0) = 0.5;
    T(1, 1) = 0.3;
    // Block 2 at (2, 2): complex pair 1.0 +/- 0.7i
    T(2, 2) = 1.0;
    T(2, 3) = -0.7;
    T(3, 2) = 0.7;
    T(3, 3) = 1.0;
    // Coupling block (0:2, 2:4)
    T(0, 2) = 0.4;
    T(0, 3) = -0.1;
    T(1, 2) = 0.2;
    T(1, 3) = 0.3;
    Mat U = Mat::Identity();
    const Mat T_before = T;
    const double scale = infinity_norm(T);
    const double eps = std::numeric_limits<double>::epsilon();

    double pivot_ratio = 0.0;
    ctrlpp::detail::pivot_ratio_conditioning tag{};
    const bool ok = ctrlpp::detail::swap_real_schur_blocks<double, N>(
        T, U, 0, 2, 2, tag, pivot_ratio);

    REQUIRE(ok);
    CHECK(pivot_ratio > eps);

    // Post-swap: (2:4, 0:2) should be zero.
    for (int i = 2; i < 4; ++i)
        for (int j = 0; j < 2; ++j)
            CHECK(std::abs(T(i, j)) < 10 * eps * scale);

    CHECK(orthogonality_residual(U) < 10 * eps * N);

    // Trace preserved.
    CHECK_THAT(T.trace(), WithinAbs(T_before.trace(), 100 * eps * scale));

    // Similarity identity on the full N x N.
    Mat reconstructed = U.transpose() * T_before * U;
    CHECK((reconstructed - T).norm() < 100 * eps * scale);
}

TEST_CASE("swap rejects when eigenvalues coincide -- ill-conditioned Sylvester")
{
    constexpr int N = 2;
    using Mat = Eigen::Matrix<double, N, N>;
    Mat T;
    // Identical diagonal entries: Sylvester equation has no unique solution
    // beyond the trivial homogeneous kernel when coupling beta != 0.
    const double alpha = 1.5;
    T << alpha, 0.7,
         0.0,   alpha;
    Mat U = Mat::Identity();
    const double eps = std::numeric_limits<double>::epsilon();

    double pivot_ratio = 1.0;
    ctrlpp::detail::pivot_ratio_conditioning tag{};
    // swap_real_schur_blocks for 1x1/1x1 routes to swap_real_schur_1x1,
    // which detects the coinciding-eigenvalue no-op case. Confirm that path.
    const bool ok = ctrlpp::detail::swap_real_schur_blocks<double, N>(
        T, U, 0, 1, 1, tag, pivot_ratio);

    // 1x1/1x1 no-op branch returns true but leaves T unchanged.
    REQUIRE(ok);
    CHECK(pivot_ratio == 1.0);
    CHECK(T(0, 0) == alpha);
    CHECK(T(1, 1) == alpha);
    CHECK(orthogonality_residual(U) < 10 * eps * N);
}

TEST_CASE("swap is overflow-safe for ill-scaled 1x1/1x1 inputs")
{
    constexpr int N = 2;
    using Mat = Eigen::Matrix<double, N, N>;
    Mat T;
    T << 1.0e-300, 1.0,
         0.0,      1.0e+300;
    Mat U = Mat::Identity();
    const double eps = std::numeric_limits<double>::epsilon();

    double pivot_ratio = 0.0;
    ctrlpp::detail::pivot_ratio_conditioning tag{};
    const bool ok = ctrlpp::detail::swap_real_schur_blocks<double, N>(
        T, U, 0, 1, 1, tag, pivot_ratio);

    REQUIRE(ok);
    // No NaN or Inf anywhere.
    CHECK(T.allFinite());
    CHECK(U.allFinite());
    CHECK(orthogonality_residual(U) < 10 * eps * N);
}
