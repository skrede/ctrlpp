#include "ctrlpp/sysid/batch_arx.h"
#include "ctrlpp/sysid.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <cstdint>
#include <algorithm>

using Catch::Matchers::WithinAbs;

namespace
{

constexpr std::size_t identity_record_length = 512;

void build_identity_record(Eigen::Matrix<double, 1, static_cast<int>(identity_record_length)>& Y,
                           Eigen::Matrix<double, 1, static_cast<int>(identity_record_length)>& U)
{
    std::uint32_t lfsr = 0x1u;
    double y1 = 0.0;
    double y2 = 0.0;
    double u1 = 0.0;
    double u2 = 0.0;
    for(std::size_t t = 0; t < identity_record_length; ++t)
    {
        std::uint32_t const feedback = ((lfsr >> 0) ^ (lfsr >> 4)) & 1u;
        lfsr = (lfsr >> 1) | (feedback << 8);
        double const u = (lfsr & 1u) ? 1.0 : -1.0;
        double const d = ((lfsr >> 3) & 1u) ? 0.03125 : -0.03125;
        double const y = 1.25 * y1 - 0.5 * y2 + 0.25 * u1 + 0.125 * u2 + d;
        Y(0, static_cast<int>(t)) = y;
        U(0, static_cast<int>(t)) = u;
        y2 = y1;
        y1 = y;
        u2 = u1;
        u1 = u;
    }
}

/// @brief Bound on how far a QR least-squares fit of this record may sit from a
/// value recorded on a different build.
///
/// Householder QR applies one reflector per parameter to the augmented
/// rows x (parameters + 1) system and then back-substitutes. A reflector touches
/// at most every entry of that array at one multiply and one add apiece, so the
/// factorization costs at most 2 * parameters * rows * (parameters + 1)
/// roundings and the back-substitution at most parameters * (parameters + 1)
/// more, each worth one unit in the last place at the scale of the value it
/// lands in. Every operation is charged whether or not it rounds and no
/// cancellation is credited, so the count bounds the departure from above rather
/// than describing it, in the manner of `ctrlpp::norm_resolution_floor`. No
/// factor here is chosen: the 2 is the multiply and the add, and the +1 is the
/// target column carried alongside the regressors.
///
/// Two builds differing only in whether Eigen's `pmadd` contracts were measured
/// four units in the last place apart on this record, so the bound clears the
/// spread it exists to absorb by about four orders of magnitude. It catches a
/// changed loop bound, a dropped term or a different factorization; it does not
/// catch a reassociation.
///
/// @cite higham2002 -- Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed., 2002, Ch. 19 (least squares)
double qr_fit_tolerance(std::size_t rows, std::size_t parameters, double scale)
{
    auto const m = static_cast<double>(rows);
    auto const n = static_cast<double>(parameters);
    auto const roundings = 2.0 * n * m * (n + 1.0) + n * (n + 1.0);
    return roundings * std::numeric_limits<double>::epsilon() * std::abs(scale);
}

}

TEST_CASE("Batch ARX reproduces its recorded fit for a fixed record", "[sysid][identity]")
{
    Eigen::Matrix<double, 1, static_cast<int>(identity_record_length)> Y;
    Eigen::Matrix<double, 1, static_cast<int>(identity_record_length)> U;
    build_identity_record(Y, U);

    auto const first = ctrlpp::batch_arx<2, 2>(Y, U);
    auto const second = ctrlpp::batch_arx<2, 2>(Y, U);
    REQUIRE(first.has_value());
    REQUIRE(second.has_value());

    // Bitwise rather than a floating matcher: two solves of one record in one
    // process share every rounding decision, so any difference at all is
    // nondeterminism and a tolerance would hide it.
    REQUIRE(first->system.A.cwiseEqual(second->system.A).all());
    REQUIRE(first->system.B.cwiseEqual(second->system.B).all());
    REQUIRE(first->metrics.nrmse == second->metrics.nrmse);
    REQUIRE(first->metrics.vaf == second->metrics.vaf);

    // The structural entries are assignments of exactly representable constants
    // and carry no rounding on any build.
    REQUIRE(first->system.A(0, 1) == 1.0);
    REQUIRE(first->system.A(1, 1) == 0.0);
    REQUIRE(first->system.C(0, 0) == 1.0);
    REQUIRE(first->system.C(0, 1) == 0.0);
    REQUIRE(first->system.D(0, 0) == 0.0);

    constexpr std::array<double, 4> recorded_theta{0x1.401c4e2ad3854p+0, -0x1.0066060368e3bp-1,
                                                   0x1.ffc684a67057dp-3, 0x1.ff99d9c085327p-4};
    constexpr double recorded_nrmse = 0x1.5ee441f237a96p-4;
    constexpr double recorded_vaf = 0x1.8d1085ace4984p+6;

    double theta_scale = 0.0;
    for(double const coefficient : recorded_theta)
        theta_scale = std::max(theta_scale, std::abs(coefficient));

    auto const regressor_rows = identity_record_length - 2;
    auto const fit_tol = qr_fit_tolerance(regressor_rows, recorded_theta.size(), theta_scale);

    REQUIRE_THAT(first->system.A(0, 0), WithinAbs(recorded_theta[0], fit_tol));
    REQUIRE_THAT(first->system.A(1, 0), WithinAbs(recorded_theta[1], fit_tol));
    REQUIRE_THAT(first->system.B(0, 0), WithinAbs(recorded_theta[2], fit_tol));
    REQUIRE_THAT(first->system.B(1, 0), WithinAbs(recorded_theta[3], fit_tol));

    // The metrics score all identity_record_length samples, not the regressor rows.
    REQUIRE_THAT(first->metrics.nrmse,
                 WithinAbs(recorded_nrmse, qr_fit_tolerance(identity_record_length, recorded_theta.size(), recorded_nrmse)));
    REQUIRE_THAT(first->metrics.vaf,
                 WithinAbs(recorded_vaf, qr_fit_tolerance(identity_record_length, recorded_theta.size(), recorded_vaf)));
}

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
    REQUIRE(result.has_value());

    REQUIRE_THAT(result->system.A(0, 0), WithinAbs(0.8, 0.01));
    REQUIRE_THAT(result->system.B(0, 0), WithinAbs(0.5, 0.01));
    REQUIRE_THAT(result->system.C(0, 0), WithinAbs(1.0, 1e-15));
    REQUIRE_THAT(result->system.D(0, 0), WithinAbs(0.0, 1e-15));

    // Fit metrics should be nearly perfect
    REQUIRE(result->metrics.nrmse < 0.01);
    REQUIRE(result->metrics.vaf > 99.0);
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
    REQUIRE(result.has_value());

    // Observer canonical form: A = [a1 1; a2 0]
    REQUIRE_THAT(result->system.A(0, 0), WithinAbs(1.2, 0.01));
    REQUIRE_THAT(result->system.A(1, 0), WithinAbs(-0.5, 0.01));

    // Fit metrics should be nearly perfect
    REQUIRE(result->metrics.nrmse < 0.01);
    REQUIRE(result->metrics.vaf > 99.0);
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
    REQUIRE(result.has_value());
    auto ss = result->system;

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
    REQUIRE(result.has_value());

    // With noise, NRMSE should be non-zero but still reasonable
    REQUIRE(result->metrics.nrmse > 0.0);
    REQUIRE(result->metrics.nrmse < 0.5);
    // VAF should be reasonable but not perfect
    REQUIRE(result->metrics.vaf > 50.0);
    REQUIRE(result->metrics.vaf < 100.0);
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
    REQUIRE(result.has_value());
    auto ss = result->system;

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
    REQUIRE(result->metrics.nrmse < 1e-9);
}
