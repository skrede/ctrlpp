#ifndef HPP_GUARD_CTRLPP_STATE_SPACE_H
#define HPP_GUARD_CTRLPP_STATE_SPACE_H

#include "ctrlpp/types.h"

#include <cstddef>

namespace ctrlpp {

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct ContinuousStateSpace {
    Matrix<Scalar, NX, NX> A;
    Matrix<Scalar, NX, NU> B;
    Matrix<Scalar, NY, NX> C;
    Matrix<Scalar, NY, NU> D;
};

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct DiscreteStateSpace {
    Matrix<Scalar, NX, NX> A;
    Matrix<Scalar, NX, NU> B;
    Matrix<Scalar, NY, NX> C;
    Matrix<Scalar, NY, NU> D;
};

template<typename S, std::size_t NX>
using SISOContinuousStateSpace = ContinuousStateSpace<S, NX, 1, 1>;

template<typename S, std::size_t NX>
using SISODiscreteStateSpace = DiscreteStateSpace<S, NX, 1, 1>;

}

#endif
