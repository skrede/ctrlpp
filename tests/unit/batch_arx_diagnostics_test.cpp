#include "ctrlpp/sysid/batch_arx.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <cstdint>

namespace
{

constexpr std::size_t record_length = 512;
constexpr std::size_t realization_order = 2;
constexpr std::size_t parameter_count = 4;
constexpr std::size_t regressor_rows = record_length - realization_order;
constexpr double disturbance_amplitude = 0.03125;

using record_row = Eigen::Matrix<double, 1, static_cast<int>(record_length)>;

/// y(t) = 1.25 y(t-1) - 0.5 y(t-2) + 0.25 u(t-1) + 0.125 u(t-2) + d(t), driven
/// by a nine-bit maximum-length shift register emitting exactly representable
/// +/-1. The target column is therefore the regressor row applied to those four
/// coefficients plus d(t), up to the recurrence's own rounding, which is what
/// lets the residual be predicted rather than observed.
void build_record(record_row& Y, record_row& U, double disturbance)
{
    std::uint32_t lfsr = 0x1u;
    double y1 = 0.0;
    double y2 = 0.0;
    double u1 = 0.0;
    double u2 = 0.0;
    for(std::size_t t = 0; t < record_length; ++t)
    {
        std::uint32_t const feedback = ((lfsr >> 0) ^ (lfsr >> 4)) & 1u;
        lfsr = (lfsr >> 1) | (feedback << 8);
        double const u = (lfsr & 1u) ? 1.0 : -1.0;
        double const d = ((lfsr >> 3) & 1u) ? disturbance : -disturbance;
        double const y = 1.25 * y1 - 0.5 * y2 + 0.25 * u1 + 0.125 * u2 + d;
        Y(0, static_cast<int>(t)) = y;
        U(0, static_cast<int>(t)) = u;
        y2 = y1;
        y1 = y;
        u2 = u1;
        u1 = u;
    }
}

/// A held input makes the two input-block columns equal entry for entry, which
/// is a rank deficiency of exactly one that no arithmetic can change. The output
/// is the shift register's own sequence rather than a model's response, because
/// the rank of the regressor depends on its columns alone and an output driven
/// by a held input settles to a constant, collapsing the output block too.
void build_held_input_record(record_row& Y, record_row& U, double level)
{
    std::uint32_t lfsr = 0x1u;
    for(std::size_t t = 0; t < record_length; ++t)
    {
        std::uint32_t const feedback = ((lfsr >> 0) ^ (lfsr >> 4)) & 1u;
        lfsr = (lfsr >> 1) | (feedback << 8);
        Y(0, static_cast<int>(t)) = (lfsr & 1u) ? 1.0 : -1.0;
        U(0, static_cast<int>(t)) = level;
    }
}

/// @brief Level below which a least-squares residual carries no information
/// about the record, only about the arithmetic that formed it.
///
/// Each residual entry is the difference of two evaluations of the same
/// `parameters`-term inner product, one accumulated when the record was
/// generated and one by the fit, at one multiply and one add per term each, so
/// at most 4 * parameters roundings at the record's scale; the Euclidean norm of
/// `rows` such entries is at most sqrt(rows) times that. Every operation is
/// charged whether or not it rounds and no cancellation is credited, and the
/// factorization's own backward error is a relative perturbation orders below
/// the margin this leaves.
///
/// @cite higham2002 -- Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed., 2002, Ch. 19 (least squares)
double residual_rounding_floor(std::size_t rows, std::size_t parameters, double scale)
{
    auto const m = static_cast<double>(rows);
    auto const n = static_cast<double>(parameters);
    return std::sqrt(m) * 4.0 * n * std::numeric_limits<double>::epsilon() * std::abs(scale);
}

/// The backward-stable rank tolerance for a Householder QR is the matrix
/// dimension times unit roundoff, applied the same way as in
/// `ctrlpp::detail::care_sign_function`. Eigen compares each pivot against this
/// multiple of the largest pivot, so the value is relative and dimensionless.
/// Eigen's own default counts the diagonal size instead, which on a tall
/// regressor sits below the rounding level the pivots actually carry.
double backward_stable_rank_threshold(std::size_t rows)
{
    return static_cast<double>(rows) * std::numeric_limits<double>::epsilon();
}

}

TEST_CASE("Batch ARX reports a full-rank fit and the row count it was formed from", "[sysid][diagnostics]")
{
    record_row Y;
    record_row U;
    build_record(Y, U, 0.0);

    auto const result = ctrlpp::batch_arx<2, 2>(Y, U);
    REQUIRE(result.has_value());

    REQUIRE(result->diagnostics.parameter_count == parameter_count);
    REQUIRE(result->diagnostics.numerical_rank == parameter_count);
    REQUIRE(result->diagnostics.effective_samples == regressor_rows);
    REQUIRE(result->diagnostics.residual_norm
            <= residual_rounding_floor(regressor_rows, parameter_count, Y.cwiseAbs().maxCoeff()));
}

TEST_CASE("Batch ARX returns a rank-deficient fit rather than refusing it", "[sysid][diagnostics]")
{
    record_row Y;
    record_row U;
    build_held_input_record(Y, U, 0.0);

    auto const unexcited = ctrlpp::batch_arx<2, 2>(Y, U);
    REQUIRE(unexcited.has_value());

    // Both input-block columns are exactly zero, so their pivots are exactly
    // zero whatever the arithmetic and whatever the threshold: the record
    // resolves its two output columns and nothing else.
    REQUIRE(unexcited->diagnostics.numerical_rank == realization_order);
    REQUIRE(unexcited->diagnostics.parameter_count == parameter_count);
    REQUIRE(unexcited->system.A.allFinite());
    REQUIRE(unexcited->system.B.allFinite());

    build_held_input_record(Y, U, 1.0);
    auto const collinear = ctrlpp::batch_arx<2, 2>(Y, U, backward_stable_rank_threshold(regressor_rows));
    REQUIRE(collinear.has_value());

    REQUIRE(collinear->diagnostics.numerical_rank == parameter_count - 1);
    REQUIRE(collinear->system.A.allFinite());
    REQUIRE(collinear->system.B.allFinite());
}

TEST_CASE("Batch ARX residual norm separates a noise-free record from a disturbed one", "[sysid][diagnostics]")
{
    record_row Y_clean;
    record_row U_clean;
    record_row Y_disturbed;
    record_row U_disturbed;
    build_record(Y_clean, U_clean, 0.0);
    build_record(Y_disturbed, U_disturbed, disturbance_amplitude);

    auto const clean = ctrlpp::batch_arx<2, 2>(Y_clean, U_clean);
    auto const disturbed = ctrlpp::batch_arx<2, 2>(Y_disturbed, U_disturbed);
    REQUIRE(clean.has_value());
    REQUIRE(disturbed.has_value());

    auto const clean_floor = residual_rounding_floor(regressor_rows, parameter_count, Y_clean.cwiseAbs().maxCoeff());
    auto const disturbed_floor =
        residual_rounding_floor(regressor_rows, parameter_count, Y_disturbed.cwiseAbs().maxCoeff());

    REQUIRE(clean->diagnostics.residual_norm <= clean_floor);
    REQUIRE(disturbed->diagnostics.residual_norm > disturbed_floor);

    // The fit returns a minimizer, so its residual cannot exceed the residual at
    // the generating coefficients, and that residual is the disturbance sequence
    // itself: sqrt(rows) times its amplitude.
    REQUIRE(disturbed->diagnostics.residual_norm
            <= std::sqrt(static_cast<double>(regressor_rows)) * disturbance_amplitude + disturbed_floor);
}

TEST_CASE("The batch ARX rank threshold changes what is reported and not what is returned", "[sysid][diagnostics]")
{
    record_row Y;
    record_row U;
    build_record(Y, U, disturbance_amplitude);

    auto const reported = ctrlpp::batch_arx<2, 2>(Y, U);
    // Eigen counts a pivot as resolved when it is STRICTLY greater than this
    // multiple of the largest pivot, so at a multiplier of one no pivot
    // qualifies, not even the largest one.
    auto const suppressed = ctrlpp::batch_arx<2, 2>(Y, U, 1.0);
    REQUIRE(reported.has_value());
    REQUIRE(suppressed.has_value());

    REQUIRE(reported->diagnostics.numerical_rank == parameter_count);
    REQUIRE(suppressed->diagnostics.numerical_rank == 0);

    // How many pivots the solve runs through is fixed during the factorization
    // from epsilon alone, so the threshold reaches the reported rank and nothing
    // else. Bitwise, because a tolerance would hide exactly the coupling this
    // assertion exists to rule out.
    REQUIRE(suppressed->system.A.cwiseEqual(reported->system.A).all());
    REQUIRE(suppressed->system.B.cwiseEqual(reported->system.B).all());
    REQUIRE(suppressed->metrics.nrmse == reported->metrics.nrmse);
    REQUIRE(suppressed->metrics.vaf == reported->metrics.vaf);
    REQUIRE(suppressed->diagnostics.residual_norm == reported->diagnostics.residual_norm);
    REQUIRE(suppressed->diagnostics.effective_samples == reported->diagnostics.effective_samples);
}
