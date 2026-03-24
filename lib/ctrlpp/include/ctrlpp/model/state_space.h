#ifndef HPP_GUARD_CTRLPP_MODEL_STATE_SPACE_H
#define HPP_GUARD_CTRLPP_MODEL_STATE_SPACE_H

#include "ctrlpp/types.h"

#include <cstddef>

namespace ctrlpp {

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct continuous_state_space
{
    Matrix<Scalar, NX, NX> A;
    Matrix<Scalar, NX, NU> B;
    Matrix<Scalar, NY, NX> C;
    Matrix<Scalar, NY, NU> D;
};

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct discrete_state_space
{
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
