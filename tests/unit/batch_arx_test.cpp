#include "ctrlpp/sysid/batch_arx.h"
#include "ctrlpp/sysid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <random>

using Catch::Matchers::WithinAbs;

TEST_CASE("Batch ARX identifies first-order SISO system")
{
    // True system: y(t) = 0.8*y(t-1) + 0.5*u(t-1)
    constexpr std::size_t N = 500;
    Eigen::Matrix<double, 1, static_cast<int>(N)> Y;
    Eigen::Matrix<double, 1, static_cast<int>(N)> U;

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> u_dist(-1.0, 1.0);

    double y = 0.0;
    double u_prev = 0.0;
    for(std::size_t t = 0; t < N; ++t)
    {
        double u = u_dist(gen);
        double y_new = 0.8 * y + 0.5 * u_prev;
        Y(0, static_cast<int>(t)) = y_new;
        U(0, static_cast<int>(t)) = u;
        y = y_new;
        u_prev = u;
    }

    auto result = ctrlpp::batch_arx<1, 1>(Y, U);

    REQUIRE_THAT(result.system.A(0, 0), WithinAbs(0.8, 0.01));
    REQUIRE_THAT(result.system.B(0, 0), WithinAbs(0.5, 0.01));
    REQUIRE_THAT(result.system.C(0, 0), WithinAbs(1.0, 1e-15));
    REQUIRE_THAT(result.system.D(0, 0), WithinAbs(0.0, 1e-15));

    // Fit metrics should be nearly perfect
    REQUIRE(result.metrics.nrmse < 0.01);
    REQUIRE(result.metrics.vaf > 99.0);
}

TEST_CASE("Batch ARX identifies second-order system")
{
    // True system: y(t) = 1.2*y(t-1) - 0.5*y(t-2) + 0.3*u(t-1) + 0.1*u(t-2)
    constexpr std::size_t N = 500;
    Eigen::Matrix<double, 1, static_cast<int>(N)> Y;
    Eigen::Matrix<double, 1, static_cast<int>(N)> U;

    std::mt19937 gen(99);
    std::uniform_real_distribution<double> u_dist(-1.0, 1.0);

    double y_prev1 = 0.0;
    double y_prev2 = 0.0;
    double u_prev1 = 0.0;
    double u_prev2 = 0.0;

    for(std::size_t t = 0; t < N; ++t)
    {
        double u = u_dist(gen);
        double y_new = 1.2 * y_prev1 - 0.5 * y_prev2 + 0.3 * u_prev1 + 0.1 * u_prev2;
        Y(0, static_cast<int>(t)) = y_new;
        U(0, static_cast<int>(t)) = u;
        y_prev2 = y_prev1;
        y_prev1 = y_new;
        u_prev2 = u_prev1;
        u_prev1 = u;
    }

    auto result = ctrlpp::batch_arx<2, 2>(Y, U);

    // Observer canonical form: A = [a1 1; a2 0]
    REQUIRE_THAT(result.system.A(0, 0), WithinAbs(1.2, 0.01));
    REQUIRE_THAT(result.system.A(1, 0), WithinAbs(-0.5, 0.01));

    // Fit metrics should be nearly perfect
    REQUIRE(result.metrics.nrmse < 0.01);
    REQUIRE(result.metrics.vaf > 99.0);
}

TEST_CASE("Batch ARX state-space simulation reproduces original data")
{
    constexpr std::size_t N = 500;
    Eigen::Matrix<double, 1, static_cast<int>(N)> Y;
    Eigen::Matrix<double, 1, static_cast<int>(N)> U;

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> u_dist(-1.0, 1.0);

    double y = 0.0;
    double u_prev = 0.0;
    for(std::size_t t = 0; t < N; ++t)
    {
        double u = u_dist(gen);
        double y_new = 0.8 * y + 0.5 * u_prev;
        Y(0, static_cast<int>(t)) = y_new;
        U(0, static_cast<int>(t)) = u;
        y = y_new;
        u_prev = u;
    }

    auto result = ctrlpp::batch_arx<1, 1>(Y, U);
    auto ss = result.system;

    // Simulate state-space model
    Eigen::Matrix<double, 1, 1> x = Eigen::Matrix<double, 1, 1>::Zero();
    double max_error = 0.0;
    for(std::size_t t = 0; t < N; ++t)
    {
        Eigen::Matrix<double, 1, 1> u_vec;
        u_vec << U(0, static_cast<int>(t));
        auto y_hat = (ss.C * x + ss.D * u_vec).eval();
        x = (ss.A * x + ss.B * u_vec).eval();
        if(t > 5)
        {
            max_error = std::max(max_error, std::abs(y_hat(0, 0) - Y(0, static_cast<int>(t))));
        }
    }
    REQUIRE(max_error < 0.05);
}

TEST_CASE("Batch ARX with noisy data produces reasonable fit")
{
    constexpr std::size_t N = 500;
    Eigen::Matrix<double, 1, static_cast<int>(N)> Y;
    Eigen::Matrix<double, 1, static_cast<int>(N)> U;

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> u_dist(-1.0, 1.0);
    std::normal_distribution<double> noise(0.0, 0.1);

    double y = 0.0;
    double u_prev = 0.0;
    for(std::size_t t = 0; t < N; ++t)
    {
        double u = u_dist(gen);
        double y_new = 0.8 * y + 0.5 * u_prev + noise(gen);
        Y(0, static_cast<int>(t)) = y_new;
        U(0, static_cast<int>(t)) = u;
        y = y_new;
        u_prev = u;
    }

    auto result = ctrlpp::batch_arx<1, 1>(Y, U);

    // With noise, NRMSE should be non-zero but still reasonable
    REQUIRE(result.metrics.nrmse > 0.0);
    REQUIRE(result.metrics.nrmse < 0.5);
    // VAF should be reasonable but not perfect
    REQUIRE(result.metrics.vaf > 50.0);
    REQUIRE(result.metrics.vaf < 100.0);
}

TEST_CASE("Batch ARX with NB > NA realizes all b-coefficients (max(NA,NB) states)")
{
    // True system: y(t) = a1*y(t-1) + b1*u(t-1) + b2*u(t-2), so NA=1 < NB=2.
    // The observer canonical realization needs max(NA,NB)=2 states; a NA-state
    // realization would drop b2 and mispredict the response from the second
    // pulse-response sample onward.
    constexpr double a1 = 0.7;
    constexpr double b1 = 0.5;
    constexpr double b2 = -0.3;

    constexpr std::size_t N = 500;
    Eigen::Matrix<double, 1, static_cast<int>(N)> Y;
    Eigen::Matrix<double, 1, static_cast<int>(N)> U;

    std::mt19937 gen(7);
    std::uniform_real_distribution<double> u_dist(-1.0, 1.0);

    double y_prev1 = 0.0;
    double u_prev1 = 0.0;
    double u_prev2 = 0.0;
    for(std::size_t t = 0; t < N; ++t)
    {
        double u = u_dist(gen);
        double y_new = a1 * y_prev1 + b1 * u_prev1 + b2 * u_prev2;
        Y(0, static_cast<int>(t)) = y_new;
        U(0, static_cast<int>(t)) = u;
        y_prev1 = y_new;
        u_prev2 = u_prev1;
        u_prev1 = u;
    }

    auto result = ctrlpp::batch_arx<1, 2>(Y, U);
    auto ss = result.system;

    // Realized state dimension is max(NA, NB) = 2, not NA = 1.
    REQUIRE(ss.A.rows() == 2);
    REQUIRE(ss.A.cols() == 2);
    REQUIRE(ss.B.rows() == 2);

    // The b2 coefficient must survive into the realization. The clearest witness
    // is the pulse response: with a unit pulse at t=0, the true ARX gives
    //   h(1) = b1,  h(2) = a1*b1 + b2,  h(3) = a1*h(2), ...
    // A realization that dropped b2 would match h(1) but diverge at h(2).
    Eigen::Matrix<double, 2, 1> x = Eigen::Matrix<double, 2, 1>::Zero();
    std::array<double, 6> h{};
    for(std::size_t k = 0; k < h.size(); ++k)
    {
        Eigen::Matrix<double, 1, 1> u_vec;
        u_vec << (k == 0 ? 1.0 : 0.0);
        h[k] = (ss.C * x + ss.D * u_vec).eval()(0, 0);
        x = (ss.A * x + ss.B * u_vec).eval();
    }

    // True pulse response of the NA=1, NB=2 plant.
    std::array<double, 6> h_true{};
    h_true[0] = 0.0;
    h_true[1] = b1;
    for(std::size_t k = 2; k < h_true.size(); ++k)
        h_true[k] = a1 * h_true[k - 1] + (k == 2 ? b2 : 0.0);

    for(std::size_t k = 0; k < h.size(); ++k)
        REQUIRE_THAT(h[k], WithinAbs(h_true[k], 1e-9));

    // The identified b2 coefficient is nonzero and recovered accurately.
    REQUIRE_THAT(ss.B(1, 0), WithinAbs(b2, 1e-9));
    REQUIRE(result.metrics.nrmse < 1e-9);
}
