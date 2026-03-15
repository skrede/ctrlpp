#ifndef HPP_GUARD_CPPCTRL_STATE_SPACE_H
#define HPP_GUARD_CPPCTRL_STATE_SPACE_H

#include "ctrlpp/linalg_policy.h"

#include <cstddef>

namespace ctrlpp {

template<LinalgPolicy Policy, typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct ContinuousStateSpace {
    typename Policy::template matrix_type<Scalar, NX, NX> A;
    typename Policy::template matrix_type<Scalar, NX, NU> B;
    typename Policy::template matrix_type<Scalar, NY, NX> C;
    typename Policy::template matrix_type<Scalar, NY, NU> D;
};

template<LinalgPolicy Policy, typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct DiscreteStateSpace {
    typename Policy::template matrix_type<Scalar, NX, NX> A;
    typename Policy::template matrix_type<Scalar, NX, NU> B;
    typename Policy::template matrix_type<Scalar, NY, NX> C;
    typename Policy::template matrix_type<Scalar, NY, NU> D;
};

template<LinalgPolicy P, typename S, std::size_t NX>
using SISOContinuousStateSpace = ContinuousStateSpace<P, S, NX, 1, 1>;

template<LinalgPolicy P, typename S, std::size_t NX>
using SISODiscreteStateSpace = DiscreteStateSpace<P, S, NX, 1, 1>;

}

#endif
