#ifndef HPP_GUARD_CTRLPP_ESTIMATION_ESTIMATION_TYPES_H
#define HPP_GUARD_CTRLPP_ESTIMATION_ESTIMATION_TYPES_H

/// @brief Public types shared by the estimation filters.
///
/// `filter_error` enumerates the structured failure modes of the fallible
/// factories in the estimation module: the quaternion attitude filters (`mekf`,
/// `manifold_ukf`, `complementary_filter`) and the sigma-point strategies and
/// the filters that build one from an options aggregate. It forms their
/// `ctrlpp::expected<T, filter_error>` construction contract.

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
enum class filter_error
{
    degenerate_quaternion,
    non_positive_sigma_spread,
    non_positive_scaling_radicand,
};

}

#endif
