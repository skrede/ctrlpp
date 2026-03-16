#ifndef HPP_GUARD_CTRLPP_PROPAGATE_H
#define HPP_GUARD_CTRLPP_PROPAGATE_H

#include "ctrlpp/state_space.h"

namespace ctrlpp {

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, LinalgPolicy Policy>
constexpr auto propagate(
    const DiscreteStateSpace<Scalar, NX, NU, NY, Policy>& sys,
    const typename Policy::template vector_type<Scalar, NX>& x,
    const typename Policy::template vector_type<Scalar, NU>& u)
    -> typename Policy::template vector_type<Scalar, NX>
{
    return Policy::add(Policy::multiply(sys.A, x), Policy::multiply(sys.B, u));
}

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, LinalgPolicy Policy>
constexpr auto output(
    const DiscreteStateSpace<Scalar, NX, NU, NY, Policy>& sys,
    const typename Policy::template vector_type<Scalar, NX>& x,
    const typename Policy::template vector_type<Scalar, NU>& u)
    -> typename Policy::template vector_type<Scalar, NY>
{
    return Policy::add(Policy::multiply(sys.C, x), Policy::multiply(sys.D, u));
}

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, LinalgPolicy Policy>
constexpr auto output(
    const ContinuousStateSpace<Scalar, NX, NU, NY, Policy>& sys,
    const typename Policy::template vector_type<Scalar, NX>& x,
    const typename Policy::template vector_type<Scalar, NU>& u)
    -> typename Policy::template vector_type<Scalar, NY>
{
    return Policy::add(Policy::multiply(sys.C, x), Policy::multiply(sys.D, u));
}

}

#endif
