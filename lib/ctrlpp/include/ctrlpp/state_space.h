#ifndef HPP_GUARD_CTRLPP_STATE_SPACE_H
#define HPP_GUARD_CTRLPP_STATE_SPACE_H

#include "ctrlpp/linalg_policy.h"

#include <cstddef>

namespace ctrlpp {

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, LinalgPolicy Policy>
struct ContinuousStateSpace {
    typename Policy::template matrix_type<Scalar, NX, NX> A;
    typename Policy::template matrix_type<Scalar, NX, NU> B;
    typename Policy::template matrix_type<Scalar, NY, NX> C;
    typename Policy::template matrix_type<Scalar, NY, NU> D;
};

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, LinalgPolicy Policy>
struct DiscreteStateSpace {
    typename Policy::template matrix_type<Scalar, NX, NX> A;
    typename Policy::template matrix_type<Scalar, NX, NU> B;
    typename Policy::template matrix_type<Scalar, NY, NX> C;
    typename Policy::template matrix_type<Scalar, NY, NU> D;
};

template<typename S, std::size_t NX, LinalgPolicy P>
using SISOContinuousStateSpace = ContinuousStateSpace<S, NX, 1, 1, P>;

template<typename S, std::size_t NX, LinalgPolicy P>
using SISODiscreteStateSpace = DiscreteStateSpace<S, NX, 1, 1, P>;

}

#endif
