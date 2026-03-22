#ifndef HPP_GUARD_CTRLPP_MPC_CONSTRAINT_MODEL_H
#define HPP_GUARD_CTRLPP_MPC_CONSTRAINT_MODEL_H

#include "ctrlpp/types.h"

#include <concepts>
#include <cstddef>

namespace ctrlpp {

/// Constrains a callable to the signature g(x, u) -> Vector<Scalar, NC>,
/// representing NC nonlinear path inequality constraints g(x,u) <= 0.
template<typename G, typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC>
concept constraint_model = requires(const G& g,
                                    const Vector<Scalar, NX>& x,
                                    const Vector<Scalar, NU>& u) {
    { g(x, u) } -> std::convertible_to<Vector<Scalar, NC>>;
};

/// Constrains a callable to the signature h(x) -> Vector<Scalar, NTC>,
/// representing NTC nonlinear terminal inequality constraints h(x_N) <= 0.
template<typename H, typename Scalar, std::size_t NX, std::size_t NTC>
concept terminal_constraint_model = requires(const H& h,
                                             const Vector<Scalar, NX>& x) {
    { h(x) } -> std::convertible_to<Vector<Scalar, NTC>>;
};

}

#endif
