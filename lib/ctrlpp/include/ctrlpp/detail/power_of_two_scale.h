#ifndef HPP_GUARD_CTRLPP_DETAIL_POWER_OF_TWO_SCALE_H
#define HPP_GUARD_CTRLPP_DETAIL_POWER_OF_TWO_SCALE_H

#include <cmath>
#include <limits>

namespace ctrlpp::detail
{

// The smallest power of two at or above a positive magnitude, clamped to the largest finite one, so dividing
// by it changes only the exponent. A zero or NaN magnitude gives one and an infinite magnitude is returned
// unchanged, so callers that decline a non-finite scale still do.
template <typename Scalar>
Scalar power_of_two_ceiling(Scalar magnitude)
{
    if(!(magnitude > Scalar{0}))
        return Scalar{1};
    if(!std::isfinite(magnitude))
        return magnitude;

    int exponent = 0;
    if(std::frexp(magnitude, &exponent) == Scalar{0.5})
        --exponent;
    if(exponent > std::numeric_limits<Scalar>::max_exponent - 1)
        exponent = std::numeric_limits<Scalar>::max_exponent - 1;
    return std::ldexp(Scalar{1}, exponent);
}

}

#endif
