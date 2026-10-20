#include "ctrlpp/estimation/estimation_types.h"
#include "ctrlpp/estimation/sigma_points/julier_sigma_points.h"
#include "ctrlpp/estimation/sigma_points/merwe_sigma_points.h"
#include "ctrlpp/estimation/sigma_points/sigma_point_strategy.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <cstddef>

using namespace ctrlpp;

// ---------------------------------------------------------------------------
// Concept satisfaction static asserts
// ---------------------------------------------------------------------------
static_assert(sigma_point_strategy<merwe_sigma_points<double, 2>, double, 2>);
static_assert(sigma_point_strategy<julier_sigma_points<double, 2>, double, 2>);

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("merwe sigma points generate correct weights")
{
    constexpr std::size_t NX = 2;
    auto sp_result = merwe_sigma_points<double, NX>::try_create(merwe_options<double>{.alpha = 1e-3, .beta = 2.0, .kappa = 0.0});
    REQUIRE(sp_result.has_value());
    auto const& sp = *sp_result;

    Vector<double, NX> x;
    x << 1.0, 2.0;
    Matrix<double, NX, NX> P = Matrix<double, NX, NX>::Identity();

    auto result = sp.generate(x, P);

    // Wm should sum to 1.0
    double wm_sum = 0.0;
    for(auto w : result.Wm)
        wm_sum += w;
    CHECK_THAT(wm_sum, Catch::Matchers::WithinAbs(1.0, 1e-10));

    // Wc sum: Wc_0 = Wm_0 + (1 - alpha^2 + beta), so Wc sums to Wm_sum + (1 - alpha^2 + beta)
    double wc_sum = 0.0;
    for(auto w : result.Wc)
        wc_sum += w;
    double expected_wc_sum = 1.0 + (1.0 - 1e-3 * 1e-3 + 2.0);
    CHECK_THAT(wc_sum, Catch::Matchers::WithinAbs(expected_wc_sum, 1e-10));

    // Verify Wm_0 and Wc_0 formulas
    double alpha = 1e-3;
    double beta = 2.0;
    double kappa = 0.0;
    double n = static_cast<double>(NX);
    double lambda = alpha * alpha * (n + kappa) - n;
    double expected_wm0 = lambda / (n + lambda);
    double expected_wc0 = expected_wm0 + (1.0 - alpha * alpha + beta);

    CHECK_THAT(result.Wm[0], Catch::Matchers::WithinAbs(expected_wm0, 1e-12));
    CHECK_THAT(result.Wc[0], Catch::Matchers::WithinAbs(expected_wc0, 1e-12));

    // Verify Wm_i = Wc_i = 1/(2*(n+lambda))
    double expected_wi = 1.0 / (2.0 * (n + lambda));
    for(std::size_t i = 1; i < merwe_sigma_points<double, NX>::num_points; ++i)
    {
        CHECK_THAT(result.Wm[i], Catch::Matchers::WithinAbs(expected_wi, 1e-12));
        CHECK_THAT(result.Wc[i], Catch::Matchers::WithinAbs(expected_wi, 1e-12));
    }
}

TEST_CASE("merwe sigma points capture mean and covariance")
{
    constexpr std::size_t NX = 2;
    auto sp_result = merwe_sigma_points<double, NX>::try_create(merwe_options<double>{.alpha = 1e-1, .beta = 2.0, .kappa = 0.0});
    REQUIRE(sp_result.has_value());
    auto const& sp = *sp_result;

    Vector<double, NX> x;
    x << 3.0, -1.0;
    Matrix<double, NX, NX> P;
    P << 4.0, 1.0, 1.0, 2.0;

    auto result = sp.generate(x, P);

    // Weighted mean should reconstruct x
    Vector<double, NX> x_recon = Vector<double, NX>::Zero();
    for(std::size_t i = 0; i < merwe_sigma_points<double, NX>::num_points; ++i)
    {
        x_recon += result.Wm[i] * result.points[i];
    }
    CHECK_THAT(x_recon(0), Catch::Matchers::WithinAbs(x(0), 1e-10));
    CHECK_THAT(x_recon(1), Catch::Matchers::WithinAbs(x(1), 1e-10));

    // Weighted covariance should reconstruct P
    Matrix<double, NX, NX> P_recon = Matrix<double, NX, NX>::Zero();
    for(std::size_t i = 0; i < merwe_sigma_points<double, NX>::num_points; ++i)
    {
        auto diff = (result.points[i] - x_recon).eval();
        P_recon += result.Wc[i] * diff * diff.transpose();
    }
    CHECK_THAT(P_recon(0, 0), Catch::Matchers::WithinAbs(P(0, 0), 1e-8));
    CHECK_THAT(P_recon(0, 1), Catch::Matchers::WithinAbs(P(0, 1), 1e-8));
    CHECK_THAT(P_recon(1, 0), Catch::Matchers::WithinAbs(P(1, 0), 1e-8));
    CHECK_THAT(P_recon(1, 1), Catch::Matchers::WithinAbs(P(1, 1), 1e-8));
}

TEST_CASE("julier sigma points satisfy concept and generate")
{
    constexpr std::size_t NX = 2;
    static_assert(sigma_point_strategy<julier_sigma_points<double, NX>, double, NX>);

    auto sp_result = julier_sigma_points<double, NX>::try_create(julier_options<double>{.kappa = 0.0});
    REQUIRE(sp_result.has_value());
    auto const& sp = *sp_result;

    Vector<double, NX> x;
    x << 1.0, 2.0;
    Matrix<double, NX, NX> P = Matrix<double, NX, NX>::Identity();

    auto result = sp.generate(x, P);

    // Wm should sum to 1.0
    double wm_sum = 0.0;
    for(auto w : result.Wm)
        wm_sum += w;
    CHECK_THAT(wm_sum, Catch::Matchers::WithinAbs(1.0, 1e-10));

    // Wc should sum to 1.0
    double wc_sum = 0.0;
    for(auto w : result.Wc)
        wc_sum += w;
    CHECK_THAT(wc_sum, Catch::Matchers::WithinAbs(1.0, 1e-10));

    // Verify Wm_0 = Wc_0 = kappa/(NX + kappa)
    double n = static_cast<double>(NX);
    double kappa = 0.0;
    double expected_w0 = kappa / (n + kappa);
    CHECK_THAT(result.Wm[0], Catch::Matchers::WithinAbs(expected_w0, 1e-12));
    CHECK_THAT(result.Wc[0], Catch::Matchers::WithinAbs(expected_w0, 1e-12));

    // Verify Wm_i = Wc_i = 1/(2*(n+kappa))
    double expected_wi = 1.0 / (2.0 * (n + kappa));
    for(std::size_t i = 1; i < julier_sigma_points<double, NX>::num_points; ++i)
    {
        CHECK_THAT(result.Wm[i], Catch::Matchers::WithinAbs(expected_wi, 1e-12));
        CHECK_THAT(result.Wc[i], Catch::Matchers::WithinAbs(expected_wi, 1e-12));
    }
}

// ---------------------------------------------------------------------------
// Merwe parameter domain
//
// The scaling term is lambda = alpha^2 (n + kappa) - n, so the weight
// denominator n + lambda and the squared sigma-point offset scale gamma^2 both
// collapse to alpha^2 (n + kappa). The divisor and the radicand are the same
// expression, which fixes the admissible domain exactly, with no tolerance:
// alpha must be finite and strictly positive, and n + kappa must be finite and
// strictly positive. Each case below asserts the specific enumerator, not
// merely that construction failed.
// ---------------------------------------------------------------------------

TEST_CASE("merwe strategy rejects a zero alpha")
{
    constexpr std::size_t NX = 2;
    auto result = merwe_sigma_points<double, NX>::try_create(merwe_options<double>{.alpha = 0.0, .beta = 2.0, .kappa = 0.0});

    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == filter_error::non_positive_sigma_spread);
}

TEST_CASE("merwe strategy rejects a negative alpha")
{
    // The case a finiteness-only oracle misses: a negative alpha yields finite
    // but wrong weights, because only alpha^2 reaches the denominator.
    constexpr std::size_t NX = 2;
    auto result = merwe_sigma_points<double, NX>::try_create(merwe_options<double>{.alpha = -1.0, .beta = 2.0, .kappa = 0.0});

    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == filter_error::non_positive_sigma_spread);
}

TEST_CASE("merwe strategy rejects a non-finite alpha")
{
    constexpr std::size_t NX = 2;

    SECTION("NaN")
    {
        auto result = merwe_sigma_points<double, NX>::try_create(
            merwe_options<double>{.alpha = std::numeric_limits<double>::quiet_NaN(), .beta = 2.0, .kappa = 0.0});

        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == filter_error::non_positive_sigma_spread);
    }

    SECTION("infinity")
    {
        auto result = merwe_sigma_points<double, NX>::try_create(
            merwe_options<double>{.alpha = std::numeric_limits<double>::infinity(), .beta = 2.0, .kappa = 0.0});

        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == filter_error::non_positive_sigma_spread);
    }
}

TEST_CASE("merwe strategy rejects a zero dimension-plus-kappa sum")
{
    constexpr std::size_t NX = 2;
    auto result = merwe_sigma_points<double, NX>::try_create(
        merwe_options<double>{.alpha = 1e-3, .beta = 2.0, .kappa = -static_cast<double>(NX)});

    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == filter_error::non_positive_scaling_radicand);
}

TEST_CASE("merwe strategy rejects a negative dimension-plus-kappa sum")
{
    constexpr std::size_t NX = 2;
    auto result = merwe_sigma_points<double, NX>::try_create(
        merwe_options<double>{.alpha = 1e-3, .beta = 2.0, .kappa = -static_cast<double>(NX) - 1.0});

    REQUIRE_FALSE(result.has_value());
    CHECK(result.error() == filter_error::non_positive_scaling_radicand);
}

TEST_CASE("merwe strategy accepts a conforming parameter set and its weights are bit-exact")
{
    constexpr std::size_t NX = 2;
    constexpr double alpha = 1e-3;
    constexpr double beta = 2.0;
    constexpr double kappa = 0.0;

    auto sp_result = merwe_sigma_points<double, NX>::try_create(merwe_options<double>{.alpha = alpha, .beta = beta, .kappa = kappa});
    REQUIRE(sp_result.has_value());

    Vector<double, NX> x;
    x << 1.0, 2.0;
    Matrix<double, NX, NX> P = Matrix<double, NX, NX>::Identity();

    auto result = sp_result->generate(x, P);

    // Same operation order as the strategy, so the comparison is exact rather
    // than tolerance-bounded: the validated construction path must not perturb
    // a single weight.
    const double n = static_cast<double>(NX);
    const double lambda = alpha * alpha * (n + kappa) - n;
    const double denom = n + lambda;

    CHECK(result.Wm[0] == lambda / denom);
    CHECK(result.Wc[0] == lambda / denom + (1.0 - alpha * alpha + beta));

    const double wi = 1.0 / (2.0 * denom);
    for(std::size_t i = 1; i < merwe_sigma_points<double, NX>::num_points; ++i)
    {
        CHECK(result.Wm[i] == wi);
        CHECK(result.Wc[i] == wi);
    }

    // The center sigma point is the mean itself, copied without arithmetic, and
    // the offsets are gamma * S columns applied to it. With P the identity, S is
    // the identity too, so each offset is exactly gamma along one axis.
    const double gamma = std::sqrt(denom);

    CHECK((result.points[0].array() == x.array()).all());
    for(std::size_t i = 0; i < NX; ++i)
    {
        Vector<double, NX> offset = Vector<double, NX>::Zero();
        offset(static_cast<int>(i)) = gamma;

        CHECK((result.points[1 + i].array() == (x + offset).array()).all());
        CHECK((result.points[1 + NX + i].array() == (x - offset).array()).all());
    }
}

// ---------------------------------------------------------------------------
// Julier parameter domain
//
// n + kappa is both the radicand of the offset scale and the weight
// denominator, so kappa must be finite and n + kappa strictly positive.
// ---------------------------------------------------------------------------

TEST_CASE("julier strategy refuses kappa outside its domain")
{
    constexpr std::size_t NX = 2;
    const double n = static_cast<double>(NX);

    auto refused = [](double kappa) {
        auto result = julier_sigma_points<double, NX>::try_create(julier_options<double>{.kappa = kappa});
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == filter_error::non_positive_scaling_radicand);
    };

    SECTION("zero dimension-plus-kappa sum") { refused(-n); }
    SECTION("negative dimension-plus-kappa sum") { refused(-n - 1.0); }
    SECTION("NaN") { refused(std::numeric_limits<double>::quiet_NaN()); }
    SECTION("infinity") { refused(std::numeric_limits<double>::infinity()); }
}

TEST_CASE("julier strategy default construction is inside the domain")
{
    constexpr std::size_t NX = 2;
    auto checked = julier_sigma_points<double, NX>::try_create();
    REQUIRE(checked.has_value());

    const Vector<double, NX> x{1.0, 2.0};
    const Matrix<double, NX, NX> P = Matrix<double, NX, NX>::Identity();
    const auto from_default = julier_sigma_points<double, NX>{}.generate(x, P);
    const auto from_factory = checked->generate(x, P);

    CHECK(from_default.Wm == from_factory.Wm);
    for(std::size_t i = 0; i < julier_sigma_points<double, NX>::num_points; ++i)
        CHECK((from_default.points[i].array() == from_factory.points[i].array()).all());
}

namespace
{

// The hex literals are recorded output of the strategy at x = (1, 2, ..., NX)
// and P = I, compared bitwise; each offset point differs from x in one
// coordinate only.
template <std::size_t NX>
void check_julier_pinned(double kappa, double w0, double wi, const std::array<double, NX>& plus, const std::array<double, NX>& minus)
{
    auto sp = julier_sigma_points<double, NX>::try_create(julier_options<double>{.kappa = kappa});
    REQUIRE(sp.has_value());

    Vector<double, NX> x;
    for(std::size_t i = 0; i < NX; ++i)
        x(static_cast<int>(i)) = static_cast<double>(i + 1);
    const auto result = sp->generate(x, Matrix<double, NX, NX>::Identity());

    CHECK((result.points[0].array() == x.array()).all());
    CHECK(result.Wm[0] == w0);
    CHECK(result.Wc[0] == w0);
    for(std::size_t i = 0; i < NX; ++i)
    {
        Vector<double, NX> up = x;
        Vector<double, NX> down = x;
        up(static_cast<int>(i)) = plus[i];
        down(static_cast<int>(i)) = minus[i];
        CHECK((result.points[1 + i].array() == up.array()).all());
        CHECK((result.points[1 + NX + i].array() == down.array()).all());
    }
    for(std::size_t i = 1; i < julier_sigma_points<double, NX>::num_points; ++i)
    {
        CHECK(result.Wm[i] == wi);
        CHECK(result.Wc[i] == wi);
    }
}

}

TEST_CASE("julier strategy output for conforming kappa is bit-exact")
{
    SECTION("kappa = 0")
    {
        check_julier_pinned<2>(0.0, 0x0p+0, 0x1p-2, {0x1.3504f333f9de6p+1, 0x1.b504f333f9de6p+1},
                               {-0x1.a827999fcef34p-2, 0x1.2bec333018866p-1});
    }

    SECTION("kappa = 1")
    {
        check_julier_pinned<2>(1.0, 0x1.5555555555555p-2, 0x1.5555555555555p-3, {0x1.5db3d742c2655p+1, 0x1.ddb3d742c2655p+1},
                               {-0x1.76cf5d0b09954p-1, 0x1.126145e9ecd58p-2});
    }

    SECTION("kappa = 3 - n with n = 4 gives a negative center weight")
    {
        check_julier_pinned<4>(-1.0, -0x1.5555555555555p-2, 0x1.5555555555555p-3,
                               {0x1.5db3d742c2655p+1, 0x1.ddb3d742c2655p+1, 0x1.2ed9eba16132ap+2, 0x1.6ed9eba16132ap+2},
                               {-0x1.76cf5d0b09954p-1, 0x1.126145e9ecd58p-2, 0x1.4498517a7b356p+0, 0x1.224c28bd3d9abp+1});
    }
}
