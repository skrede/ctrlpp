#ifndef HPP_GUARD_CTRLPP_MPC_DYNAMICS_MODEL_H
#define HPP_GUARD_CTRLPP_MPC_DYNAMICS_MODEL_H

#include "ctrlpp/types.h"

#include <concepts>
#include <cstddef>

namespace ctrlpp {

/// Constrains a callable to the signature f(x, u) -> x_next,
/// where x is Vector<Scalar, NX> and u is Vector<Scalar, NU>.
template<typename D, typename Scalar, std::size_t NX, std::size_t NU>
concept dynamics_model = requires(const D& d,
                                  const Vector<Scalar, NX>& x,
                                  const Vector<Scalar, NU>& u) {
    { d(x, u) } -> std::convertible_to<Vector<Scalar, NX>>;
};

}

#endif
