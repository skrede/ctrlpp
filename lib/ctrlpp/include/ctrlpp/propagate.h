#ifndef HPP_GUARD_CPPCTRL_PROPAGATE_H
#define HPP_GUARD_CPPCTRL_PROPAGATE_H

#include "ctrlpp/state_space.h"

namespace ctrlpp {

template<LinalgPolicy Policy, typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
constexpr auto propagate(
    const DiscreteStateSpace<Policy, Scalar, NX, NU, NY>& sys,
    const typename Policy::template vector_type<Scalar, NX>& x,
    const typename Policy::template vector_type<Scalar, NU>& u)
    -> typename Policy::template vector_type<Scalar, NX>
{
    return Policy::add(Policy::multiply(sys.A, x), Policy::multiply(sys.B, u));
}

template<LinalgPolicy Policy, typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
constexpr auto output(
    const DiscreteStateSpace<Policy, Scalar, NX, NU, NY>& sys,
    const typename Policy::template vector_type<Scalar, NX>& x,
    const typename Policy::template vector_type<Scalar, NU>& u)
    -> typename Policy::template vector_type<Scalar, NY>
{
    return Policy::add(Policy::multiply(sys.C, x), Policy::multiply(sys.D, u));
}

template<LinalgPolicy Policy, typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
constexpr auto output(
    const ContinuousStateSpace<Policy, Scalar, NX, NU, NY>& sys,
    const typename Policy::template vector_type<Scalar, NX>& x,
    const typename Policy::template vector_type<Scalar, NU>& u)
    -> typename Policy::template vector_type<Scalar, NY>
{
    return Policy::add(Policy::multiply(sys.C, x), Policy::multiply(sys.D, u));
}

}

#endif
