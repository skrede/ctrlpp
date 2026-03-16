#ifndef HPP_GUARD_CTRLPP_TRANSFER_FUNCTION_H
#define HPP_GUARD_CTRLPP_TRANSFER_FUNCTION_H

#include "ctrlpp/linalg_policy.h"

#include <array>
#include <cstddef>

namespace ctrlpp {

template<typename Scalar, std::size_t NumDeg, std::size_t DenDeg, LinalgPolicy Policy>
struct TransferFunction {
    std::array<Scalar, NumDeg + 1> numerator;
    std::array<Scalar, DenDeg + 1> denominator;
};

}

#endif
