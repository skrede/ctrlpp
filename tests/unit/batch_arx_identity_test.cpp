#include "ctrlpp/sysid/batch_arx.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <limits>
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
    REQUIRE(first->diagnostics.residual_norm == second->diagnostics.residual_norm);

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
