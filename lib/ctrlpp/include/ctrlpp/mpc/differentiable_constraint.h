#ifndef HPP_GUARD_CTRLPP_MPC_DIFFERENTIABLE_CONSTRAINT_H
#define HPP_GUARD_CTRLPP_MPC_DIFFERENTIABLE_CONSTRAINT_H

#include "ctrlpp/model/constraint_model.h"

#include <cstddef>

namespace ctrlpp
{

/// Refines constraint_model by additionally requiring an analytic Jacobian
/// dg/d[x;u] returning a combined NC x (NX+NU) matrix.
template <typename G, typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC>
concept differentiable_constraint = constraint_model<G, Scalar, NX, NU, NC> && requires(const G& g, const Vector<Scalar, NX>& x, const Vector<Scalar, NU>& u) {
    { g.jacobian(x, u) } -> std::convertible_to<Matrix<Scalar, NC, NX + NU>>;
};

/// Refines terminal_constraint_model by additionally requiring an analytic
/// Jacobian dh/dx returning an NTC x NX matrix.
template <typename H, typename Scalar, std::size_t NX, std::size_t NTC>
concept differentiable_terminal_constraint = terminal_constraint_model<H, Scalar, NX, NTC> && requires(const H& h, const Vector<Scalar, NX>& x) {
    { h.jacobian(x) } -> std::convertible_to<Matrix<Scalar, NTC, NX>>;
};

} // namespace ctrlpp

#endif
