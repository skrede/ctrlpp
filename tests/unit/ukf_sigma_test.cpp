#include "ctrlpp/estimation/sigma_points/julier_sigma_points.h"
#include "ctrlpp/estimation/sigma_points/merwe_sigma_points.h"
#include "ctrlpp/estimation/sigma_points/sigma_point_strategy.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

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
    merwe_sigma_points<double, NX> sp{merwe_options<double>{.alpha = 1e-3, .beta = 2.0, .kappa = 0.0}};

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
    merwe_sigma_points<double, NX> sp{merwe_options<double>{.alpha = 1e-1, .beta = 2.0, .kappa = 0.0}};

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

    julier_sigma_points<double, NX> sp{julier_options<double>{.kappa = 0.0}};

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
