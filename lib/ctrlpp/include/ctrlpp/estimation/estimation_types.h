#ifndef HPP_GUARD_CTRLPP_ESTIMATION_ESTIMATION_TYPES_H
#define HPP_GUARD_CTRLPP_ESTIMATION_ESTIMATION_TYPES_H

/// @brief Public types shared by the estimation filters.
///
/// `filter_error` enumerates the structured failure modes of the fallible
/// factories on the quaternion attitude filters (`mekf`, `manifold_ukf`,
/// `complementary_filter`), forming their
/// `ctrlpp::expected<Filter, filter_error>` construction contract.

namespace ctrlpp
{

/// @brief Structured failure modes for the quaternion filter factories.
///
///  * degenerate_quaternion : the initial quaternion q0 has zero or non-finite
///                            norm, so normalizing it would produce NaN and
///                            silently poison the filter state at construction.
///                            A rotation quaternion is nonzero by definition;
///                            any finite nonzero q0 normalizes to a valid unit
///                            quaternion and is accepted.
enum class filter_error
{
    degenerate_quaternion,
};

}

#endif
