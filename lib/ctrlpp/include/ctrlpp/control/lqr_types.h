#ifndef HPP_GUARD_CTRLPP_CONTROL_LQR_TYPES_H
#define HPP_GUARD_CTRLPP_CONTROL_LQR_TYPES_H

#include "ctrlpp/util/concepts.h"

#include <Eigen/Core>

#include <cstddef>

namespace ctrlpp
{

// Every enumerator is an exact precondition of indexing or of the arithmetic's domain, so none carries a tolerance.
enum class lqr_error
{
    sequence_length_mismatch,   // the four lqr_tv_gains sequences hold different counts
    trajectory_length_mismatch, // lqr_cost states number neither the inputs nor one more
    step_out_of_range,          // the lqr_time_varying step is at or past the horizon
    non_finite_input,           // a matrix or vector entry is NaN or infinite
};

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct lqi_result
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NY > 0, "Output dimension NY must be positive");

    Eigen::Matrix<Scalar, int(NU), int(NX)> Kx;
    Eigen::Matrix<Scalar, int(NU), int(NY)> Ki;
};

}

#endif
