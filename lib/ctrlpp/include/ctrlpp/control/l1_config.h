#ifndef HPP_GUARD_CTRLPP_CONTROL_L1_CONFIG_H
#define HPP_GUARD_CTRLPP_CONTROL_L1_CONFIG_H

#include "ctrlpp/types.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/model/state_space.h"

#include <limits>
#include <cstddef>

namespace ctrlpp
{

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU>
struct l1_config
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    discrete_state_space<Scalar, NX, NU, NX> predictor_model{};
    Matrix<Scalar, NU, NU> gamma = Matrix<Scalar, NU, NU>::Identity();
    // Projection bounds default to an unbounded range so a default-constructed
    // config adapts freely. Set finite bounds to enable elementwise projection.
    Vector<Scalar, NU> theta_min =
        Vector<Scalar, NU>::Constant(-std::numeric_limits<Scalar>::infinity());
    Vector<Scalar, NU> theta_max =
        Vector<Scalar, NU>::Constant(std::numeric_limits<Scalar>::infinity());
    Vector<Scalar, NX> x_hat_0 = Vector<Scalar, NX>::Zero();
    Vector<Scalar, NU> sigma_hat_0 = Vector<Scalar, NU>::Zero();
};

}

#endif
