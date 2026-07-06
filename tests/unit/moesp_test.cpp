#include "ctrlpp/sysid/moesp.h"
#include "ctrlpp/sysid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>

#include <cmath>
#include <random>

using Catch::Matchers::WithinAbs;

namespace
{

// Known 2nd-order discrete state-space system (discretized spring-mass-damper)
// Eigenvalues inside unit circle for stability
constexpr double a11 = 0.9;
constexpr double a12 = 0.1;
constexpr double a21 = -0.1;
constexpr double a22 = 0.85;

struct test_data
{
    Eigen::Matrix<double, 1, Eigen::Dynamic> Y;
    Eigen::Matrix<double, 1, Eigen::Dynamic> U;
};

auto generate_data(std::size_t N, double noise_std = 0.0) -> test_data
{
    Eigen::Matrix2d A;
    A << a11, a12, a21, a22;
    Eigen::Vector2d B;
    B << 0.0, 0.1;
    Eigen::RowVector2d C;
    C << 1.0, 0.0;

    Eigen::Matrix<double, 1, Eigen::Dynamic> Y(1, static_cast<Eigen::Index>(N));
    Eigen::Matrix<double, 1, Eigen::Dynamic> U(1, static_cast<Eigen::Index>(N));

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> u_dist(-1.0, 1.0);
    Eigen::Vector2d x = Eigen::Vector2d::Zero();
    for(std::size_t t = 0; t < N; ++t)
    {
        double u = u_dist(gen);
        U(0, static_cast<Eigen::Index>(t)) = u;

        double y = (C * x)(0);
        if(noise_std > 0.0)
        {
            std::normal_distribution<double> noise(0.0, noise_std);
            y += noise(gen);
        }
        Y(0, static_cast<Eigen::Index>(t)) = y;

        x = A * x + B * u;
    }

    return {Y, U};
}

} // namespace

TEST_CASE("MOESP singular values show clear gap for 2nd-order system")
{
    auto [Y, U] = generate_data(1500);

    auto sv = ctrlpp::moesp_singular_values(Y, U);

    REQUIRE(sv.size() >= 2);
    // First two should be significantly larger than the rest
    REQUIRE(sv(0) > sv(1));
    // Gap: ratio of 2nd to 3rd singular value should be large
    if(sv.size() > 2)
    {
        double ratio = sv(1) / sv(2);
        REQUIRE(ratio > 5.0);
    }
}

TEST_CASE("MOESP identifies 2nd-order system with good fit")
{
    auto [Y, U] = generate_data(1500);

    auto result = ctrlpp::moesp<2>(Y, U);

    // Simulate identified model on input data and check output match
    REQUIRE(result.metrics.nrmse < 0.1);
}

TEST_CASE("MOESP condition number is finite and positive")
{
    auto [Y, U] = generate_data(1500);

    auto result = ctrlpp::moesp<2>(Y, U);

    REQUIRE(result.condition_number > 0.0);
    REQUIRE(std::isfinite(result.condition_number));
}

TEST_CASE("MOESP VAF > 90 for clean data")
{
    auto [Y, U] = generate_data(1500);

    auto result = ctrlpp::moesp<2>(Y, U);

    REQUIRE(result.metrics.vaf > 90.0);
}

TEST_CASE("MOESP with measurement noise produces degraded but positive VAF")
{
    auto [Y, U] = generate_data(1500, 0.05);

    auto result = ctrlpp::moesp<2>(Y, U);

    REQUIRE(result.metrics.vaf > 0.0);
    // Should still identify something reasonable with mild noise
    REQUIRE(result.metrics.vaf < 100.0);
}

// The output row C must be the first block row of the observability matrix
// Gamma, not its first column (which only coincides for NX = 1). A wrong C is
// invisible to fit metrics, eig(A), and Markov parameters, because the B/D
// least-squares refit produces a valid similar realization that reproduces the
// transfer function exactly. It is caught only by the subspace geometry: for the
// correct realization the reconstructed observability matrix O = [C; CA; ...] is
// exactly Gamma = U * sqrt(S), whose Gram is O^T O = diag(singular_values).
TEST_CASE("MOESP recovers C consistent with the observability subspace")
{
    auto [Y, U] = generate_data(1500);

    constexpr int i_blocks = 10;
    auto result = ctrlpp::moesp<2>(Y, U, i_blocks);
    auto const& sys = result.system;

    // The identified poles match the true plant: the characteristic polynomial
    // (trace, determinant) is similarity invariant. True A has trace 1.75 and
    // determinant 0.775.
    Eigen::Matrix2d A_id = sys.A;
    REQUIRE_THAT(A_id.trace(), WithinAbs(a11 + a22, 1e-6));
    REQUIRE_THAT(A_id.determinant(), WithinAbs(a11 * a22 - a12 * a21, 1e-6));

    // Reconstruct the observability matrix O = [C; CA; ...; CA^{i-1}] from the
    // returned (A, C). For the correct C it equals Gamma, so O^T O = diag(sv);
    // extracting C from Gamma's first column instead breaks this (O^T O gains
    // large off-diagonal terms and its diagonal no longer matches sv).
    Eigen::MatrixXd O(i_blocks, 2);
    Eigen::RowVector2d row = sys.C;
    for(int k = 0; k < i_blocks; ++k)
    {
        O.row(k) = row;
        row = row * A_id;
    }
    Eigen::Matrix2d OtO = O.transpose() * O;

    REQUIRE(result.singular_values.size() >= 2);
    double const sv0 = result.singular_values(0);
    double const sv1 = result.singular_values(1);
    double const scale = sv0; // largest singular value sets the magnitude

    REQUIRE_THAT(OtO(0, 0), WithinAbs(sv0, 1e-9 * scale));
    REQUIRE_THAT(OtO(1, 1), WithinAbs(sv1, 1e-9 * scale));
    REQUIRE_THAT(OtO(0, 1), WithinAbs(0.0, 1e-9 * scale));
    REQUIRE_THAT(OtO(1, 0), WithinAbs(0.0, 1e-9 * scale));

    // The realization also reproduces the true impulse response (Markov params).
    Eigen::Matrix2d A_true;
    A_true << a11, a12, a21, a22;
    Eigen::Vector2d B_true(0.0, 0.1);
    Eigen::RowVector2d C_true(1.0, 0.0);

    Eigen::Vector2d x_true = Eigen::Vector2d::Zero();
    Eigen::Vector2d x_id = Eigen::Vector2d::Zero();
    for(int k = 0; k < 8; ++k)
    {
        double const u = (k == 0) ? 1.0 : 0.0;
        double const h_true = (C_true * x_true)(0);
        double const h_id = (sys.C * x_id)(0) + sys.D(0, 0) * u;
        REQUIRE_THAT(h_id, WithinAbs(h_true, 1e-9));
        x_true = A_true * x_true + B_true * u;
        Eigen::Matrix<double, 1, 1> uv;
        uv(0, 0) = u;
        x_id = (sys.A * x_id + sys.B * uv).eval();
    }
}
