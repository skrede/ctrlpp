#ifndef HPP_GUARD_CTRLPP_MODEL_DIFFERENTIABLE_DYNAMICS_H
#define HPP_GUARD_CTRLPP_MODEL_DIFFERENTIABLE_DYNAMICS_H

#include "ctrlpp/model/dynamics_model.h"

#include <cstddef>

namespace ctrlpp
{

/// Refines dynamics_model by additionally requiring analytic Jacobians
/// with respect to state (jacobian_x) and input (jacobian_u).
template <typename D, typename Scalar, std::size_t NX, std::size_t NU>
concept differentiable_dynamics = dynamics_model<D, Scalar, NX, NU> && requires(const D& d, const Vector<Scalar, NX>& x, const Vector<Scalar, NU>& u) {
    { d.jacobian_x(x, u) } -> std::convertible_to<Matrix<Scalar, NX, NX>>;
    { d.jacobian_u(x, u) } -> std::convertible_to<Matrix<Scalar, NX, NU>>;
};

} // namespace ctrlpp

#endif
