#ifndef HPP_GUARD_CTRLPP_MODEL_MEASUREMENT_MODEL_H
#define HPP_GUARD_CTRLPP_MODEL_MEASUREMENT_MODEL_H

#include "ctrlpp/types.h"

#include <concepts>
#include <cstddef>

namespace ctrlpp {

/// Constrains a callable to the signature h(x) -> y,
/// where x is Vector<Scalar, NX> and y is Vector<Scalar, NY>.
template<typename M, typename Scalar, std::size_t NX, std::size_t NY>
concept measurement_model = requires(const M& m,
                                     const Vector<Scalar, NX>& x) {
    { m(x) } -> std::convertible_to<Vector<Scalar, NY>>;
};

}

#endif
