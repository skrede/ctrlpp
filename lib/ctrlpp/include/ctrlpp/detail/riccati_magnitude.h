#ifndef HPP_GUARD_CTRLPP_DETAIL_RICCATI_MAGNITUDE_H
#define HPP_GUARD_CTRLPP_DETAIL_RICCATI_MAGNITUDE_H

#include <Eigen/Core>

#include <cmath>
#include <limits>
#include <cstdint>
#include <initializer_list>

namespace ctrlpp::detail
{

// A magnitude whose `value` means something only when `resolved`. Unresolved is contagious and no comparison with an
// unresolved side succeeds, because a plain Frobenius norm overflows to infinity (where inf <= inf certifies anything)
// or underflows to zero (where 0 <= 0 passes): a magnitude the arithmetic could not form declines instead.
template <typename Scalar>
struct resolved_magnitude
{
    Scalar value;
    bool resolved;
};

// The largest-entry rescale max|M| * ||M / max|M|||, which cannot leave the range before its final product. Eigen's
// blueNorm() returns zero for the smallest subnormal, and stableNorm() costs more at the sizes both Riccati paths carry.
template <typename Derived>
resolved_magnitude<typename Derived::Scalar> resolve_magnitude(const Eigen::MatrixBase<Derived>& operand)
{
    using Scalar = typename Derived::Scalar;
    const auto& evaluated = operand.eval();

    if(!evaluated.allFinite())
        return {Scalar{0}, false};

    const Scalar largest_entry = evaluated.cwiseAbs().maxCoeff();
    if(!(largest_entry > Scalar{0}))
        return {Scalar{0}, true};

    const Scalar magnitude = largest_entry * (evaluated / largest_entry).norm();
    if(!std::isfinite(magnitude))
        return {Scalar{0}, false};
    return {magnitude, true};
}

// Squares nothing, so only the finiteness guard is needed. A bound derived from the largest entry differs from the
// Frobenius one, and substituting one for the other would move a threshold.
template <typename Derived>
resolved_magnitude<typename Derived::Scalar> resolve_largest_entry(const Eigen::MatrixBase<Derived>& operand)
{
    using Scalar = typename Derived::Scalar;
    const auto& evaluated = operand.eval();

    if(!evaluated.allFinite())
        return {Scalar{0}, false};
    return {evaluated.cwiseAbs().maxCoeff(), true};
}

template <typename Scalar>
resolved_magnitude<Scalar> largest_magnitude(std::initializer_list<resolved_magnitude<Scalar>> operands)
{
    resolved_magnitude<Scalar> largest{Scalar{0}, true};
    for(const resolved_magnitude<Scalar>& operand : operands)
    {
        if(!operand.resolved)
            return {Scalar{0}, false};
        if(operand.value > largest.value)
            largest.value = operand.value;
    }
    return largest;
}

template <typename Scalar>
resolved_magnitude<Scalar> scaled_magnitude(const resolved_magnitude<Scalar>& operand, Scalar factor)
{
    if(!operand.resolved || !std::isfinite(factor))
        return {Scalar{0}, false};

    const Scalar scaled = operand.value * factor;
    if(!std::isfinite(scaled))
        return {Scalar{0}, false};
    return {scaled, true};
}

template <typename Scalar>
resolved_magnitude<Scalar> counted_rounding_bound(const resolved_magnitude<Scalar>& scale, std::int32_t rounding_ops)
{
    return scaled_magnitude(scale, static_cast<Scalar>(rounding_ops) * std::numeric_limits<Scalar>::epsilon());
}

template <typename Scalar>
bool magnitude_within(const resolved_magnitude<Scalar>& value, const resolved_magnitude<Scalar>& bound)
{
    return value.resolved && bound.resolved && value.value <= bound.value;
}

}

#endif
