#ifndef HPP_GUARD_CTRLPP_CONTROL_L1_CONFIG_H
#define HPP_GUARD_CTRLPP_CONTROL_L1_CONFIG_H

#include "ctrlpp/types.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/model/state_space.h"

#include <limits>
#include <cstddef>

namespace ctrlpp
{

enum class l1_error
{
    singular_predictor,
    singular_dc_gain,
    non_finite_gain,
    invalid_filter_config,
};

// Listed in the order l1_controller::evaluate tests them, which follows the controller's data flow, so the cause named
// is the most upstream one.
enum class l1_step_error
{
    non_finite_predictor_state,
    non_finite_adaptation,
    non_finite_state,
    non_finite_reference,
    non_finite_result,
};

// Declared in severity order: the controller only ever raises its status, so a lesser fault cannot mask a greater one.
enum class l1_health
{
    ok,
    projection_clamped_non_finite,
    non_finite_carried_state,
};

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU>
struct l1_config
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    discrete_state_space<Scalar, NX, NU, NX> predictor_model{};
    Matrix<Scalar, NU, NU> gamma = Matrix<Scalar, NU, NU>::Identity();
    Vector<Scalar, NU> theta_min = Vector<Scalar, NU>::Constant(-std::numeric_limits<Scalar>::infinity());
    Vector<Scalar, NU> theta_max = Vector<Scalar, NU>::Constant(std::numeric_limits<Scalar>::infinity());
    Vector<Scalar, NX> x_hat_0 = Vector<Scalar, NX>::Zero();
    Vector<Scalar, NU> sigma_hat_0 = Vector<Scalar, NU>::Zero();
};

}

#endif
