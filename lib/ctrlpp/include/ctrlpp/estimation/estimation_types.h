#ifndef HPP_GUARD_CTRLPP_ESTIMATION_ESTIMATION_TYPES_H
#define HPP_GUARD_CTRLPP_ESTIMATION_ESTIMATION_TYPES_H

/// @brief Public types shared by the estimation filters.
///
/// `filter_error` enumerates the structured failure modes of the fallible
/// factories in the estimation module: the quaternion attitude filters (`mekf`,
/// `manifold_ukf`, `complementary_filter`), the sigma-point strategies, the
/// filters that build one from an options aggregate, and the three Kalman
/// filters' configuration validation. It forms their
/// `ctrlpp::expected<T, filter_error>` construction contract.

#include "ctrlpp/expected.h"

namespace ctrlpp
{

/// @brief Structured failure modes for the estimation-module factories.
///
///  * degenerate_quaternion         : the initial quaternion q0 has zero or
///                                    non-finite norm, so normalizing it would
///                                    produce NaN and silently poison the filter
///                                    state at construction. A rotation
///                                    quaternion is nonzero by definition; any
///                                    finite nonzero q0 normalizes to a valid
///                                    unit quaternion and is accepted.
///  * non_positive_sigma_spread     : the sigma-point spread parameter is not
///                                    finite and strictly positive. It is the
///                                    divisor of the sigma-point weight
///                                    denominator, so this is an exact domain
///                                    condition and not a tuning preference: a
///                                    zero spread makes the weights infinite,
///                                    and a negative one makes them finite but
///                                    wrong, since only its square is used.
///  * non_positive_scaling_radicand : the sum of the state dimension and the
///                                    secondary scaling parameter is not finite
///                                    and strictly positive. It sits under the
///                                    square root that scales the sigma-point
///                                    offsets and inside the same weight
///                                    denominator, so a non-positive sum yields
///                                    non-finite offsets or non-finite weights.
///  * non_finite_process_noise      : the process-noise matrix Q has a
///                                    non-finite entry. Q is added to the
///                                    propagated covariance every prediction, so
///                                    an infinite entry makes the predicted
///                                    covariance infinite, the gain a ratio of
///                                    infinities and therefore non-finite, and
///                                    the corrected estimate non-finite with it.
///  * non_finite_measurement_noise  : the measurement-noise matrix R has a
///                                    non-finite entry. R is added to the
///                                    innovation covariance, which is the matrix
///                                    the gain solve is posed against, so a
///                                    non-finite entry makes that solve
///                                    meaningless for every measurement.
///  * non_finite_initial_state      : the initial state vector x0 has a
///                                    non-finite component. It seeds the carried
///                                    estimate, which is the filter's memory, so
///                                    every later estimate is formed from it.
///  * non_finite_initial_covariance : the initial covariance matrix P0 has a
///                                    non-finite entry. It seeds the covariance
///                                    recursion, and the recursion has no
///                                    mechanism that could return a non-finite
///                                    covariance to a finite one.
///
/// The four configuration enumerators state finiteness and nothing else. An
/// ill-conditioned but finite configuration -- a covariance with entries many
/// orders of magnitude apart, or a singular P0 -- is a legitimate, deliberately
/// posed problem and is accepted: conditioning is a numerical-behavior question
/// and finiteness is the domain condition. Symmetry and positive definiteness
/// are likewise not tested here, because the filters do not promise them of the
/// configuration they are handed.
enum class filter_error
{
    degenerate_quaternion,
    non_positive_sigma_spread,
    non_positive_scaling_radicand,
    non_finite_process_noise,
    non_finite_measurement_noise,
    non_finite_initial_state,
    non_finite_initial_covariance,
};

namespace detail
{

/// @brief Validate the noise and initial-condition fields that `kalman_config`,
/// `ekf_config` and `ukf_config` declare identically.
///
/// The three aggregates carry the same four fields in the same order and feed
/// them into the same recursion, so the condition is shared rather than
/// per-filter and lives here rather than being copied three times.
///
/// The fields are independent of one another -- unlike a step's carried estimate
/// and its measurement, no one of them is upstream of another -- so the report
/// order is the order the aggregates declare, not a claim about causality.
///
/// Finiteness is the whole test. Nothing here rejects a small, singular or
/// asymmetric matrix.
template <typename Process, typename Measurement, typename State, typename Covariance>
auto validate_filter_configuration(const Process& Q, const Measurement& R, const State& x0, const Covariance& P0) -> ctrlpp::expected<void, filter_error>
{
    if(!Q.allFinite())
        return ctrlpp::unexpected(filter_error::non_finite_process_noise);
    if(!R.allFinite())
        return ctrlpp::unexpected(filter_error::non_finite_measurement_noise);
    if(!x0.allFinite())
        return ctrlpp::unexpected(filter_error::non_finite_initial_state);
    if(!P0.allFinite())
        return ctrlpp::unexpected(filter_error::non_finite_initial_covariance);
    return {};
}

}

}

#endif
