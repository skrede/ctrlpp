#ifndef HPP_GUARD_CTRLPP_MODEL_DIFFERENTIABLE_MEASUREMENT_H
#define HPP_GUARD_CTRLPP_MODEL_DIFFERENTIABLE_MEASUREMENT_H

#include "ctrlpp/model/measurement_model.h"

#include <cstddef>

namespace ctrlpp
{

/// Refines measurement_model by additionally requiring an analytic Jacobian
/// of the measurement function with respect to state.
template <typename M, typename Scalar, std::size_t NX, std::size_t NY>
concept differentiable_measurement = measurement_model<M, Scalar, NX, NY> && requires(const M& m, const Vector<Scalar, NX>& x) {
    { m.jacobian(x) } -> std::convertible_to<Matrix<Scalar, NY, NX>>;
};

} // namespace ctrlpp

#endif
