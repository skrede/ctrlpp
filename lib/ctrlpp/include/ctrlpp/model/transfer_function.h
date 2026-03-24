#ifndef HPP_GUARD_CTRLPP_MODEL_TRANSFER_FUNCTION_H
#define HPP_GUARD_CTRLPP_MODEL_TRANSFER_FUNCTION_H

/// @brief SISO transfer function representation (numerator/denominator polynomial).
///
/// @cite franklin2015 -- Franklin et al., "Feedback Control of Dynamic Systems", 2015

#include <array>
#include <cstddef>

namespace ctrlpp
{

template <typename Scalar, std::size_t NumDeg, std::size_t DenDeg>
struct transfer_function
{
    std::array<Scalar, NumDeg + 1> numerator;
    std::array<Scalar, DenDeg + 1> denominator;
};

} // namespace ctrlpp

#endif
