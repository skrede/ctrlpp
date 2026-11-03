#ifndef HPP_GUARD_CTRLPP_ESTIMATION_UKF_TYPES_H
#define HPP_GUARD_CTRLPP_ESTIMATION_UKF_TYPES_H

#include "ctrlpp/types.h"

#include "ctrlpp/util/concepts.h"

#include <cstddef>

namespace ctrlpp
{

enum class gain_decomposition
{
    ldlt,
    qr
};

enum class ukf_update_error
{
    non_finite_state,
    non_finite_covariance,
    non_finite_measurement,
    non_finite_result,
};

enum class ukf_health
{
    ok,
    covariance_repaired,
    non_finite_estimate
};

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct ukf_config
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NY > 0, "Output dimension NY must be positive");
    Matrix<Scalar, NX, NX> Q{Matrix<Scalar, NX, NX>::Identity()};
    Matrix<Scalar, NY, NY> R{Matrix<Scalar, NY, NY>::Identity()};
    Vector<Scalar, NX> x0{Vector<Scalar, NX>::Zero()};
    Matrix<Scalar, NX, NX> P0{Matrix<Scalar, NX, NX>::Identity()};
    gain_decomposition decomposition{gain_decomposition::ldlt};
};

}

#endif
