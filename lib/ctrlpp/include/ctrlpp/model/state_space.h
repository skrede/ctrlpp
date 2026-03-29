#ifndef HPP_GUARD_CTRLPP_MODEL_STATE_SPACE_H
#define HPP_GUARD_CTRLPP_MODEL_STATE_SPACE_H

#include "ctrlpp/types.h"

#include <cstddef>
#include <type_traits>

namespace ctrlpp
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct continuous_state_space
{
    static_assert(std::is_floating_point_v<Scalar>,
                  "continuous_state_space requires a floating-point Scalar type");
    static_assert(NX >= 1 && NU >= 1 && NY >= 1,
                  "state_space dimensions must be >= 1");

    Matrix<Scalar, NX, NX> A;
    Matrix<Scalar, NX, NU> B;
    Matrix<Scalar, NY, NX> C;
    Matrix<Scalar, NY, NU> D;
};

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct discrete_state_space
{
    static_assert(std::is_floating_point_v<Scalar>,
                  "discrete_state_space requires a floating-point Scalar type");
    static_assert(NX >= 1 && NU >= 1 && NY >= 1,
                  "state_space dimensions must be >= 1");

    Matrix<Scalar, NX, NX> A;
    Matrix<Scalar, NX, NU> B;
    Matrix<Scalar, NY, NX> C;
    Matrix<Scalar, NY, NU> D;
};

template <typename S, std::size_t NX>
using siso_continuous_state_space = continuous_state_space<S, NX, 1, 1>;

template <typename S, std::size_t NX>
using siso_discrete_state_space = discrete_state_space<S, NX, 1, 1>;

}

#endif
