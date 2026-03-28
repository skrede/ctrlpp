#ifndef HPP_GUARD_CTRLPP_CONTROL_L1_CONFIG_H
#define HPP_GUARD_CTRLPP_CONTROL_L1_CONFIG_H

#include "ctrlpp/types.h"

#include "ctrlpp/model/state_space.h"

#include <cstddef>

namespace ctrlpp
{

template <typename Scalar, std::size_t NX, std::size_t NU>
struct l1_config
{
    discrete_state_space<Scalar, NX, NU, NX> predictor_model{};
    Matrix<Scalar, NU, NU> gamma = Matrix<Scalar, NU, NU>::Identity();
    Vector<Scalar, NU> theta_min = Vector<Scalar, NU>::Zero();
    Vector<Scalar, NU> theta_max = Vector<Scalar, NU>::Zero();
    Vector<Scalar, NX> x_hat_0 = Vector<Scalar, NX>::Zero();
    Vector<Scalar, NU> sigma_hat_0 = Vector<Scalar, NU>::Zero();
};

}

#endif
