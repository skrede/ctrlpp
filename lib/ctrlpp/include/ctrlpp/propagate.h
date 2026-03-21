#ifndef HPP_GUARD_CTRLPP_PROPAGATE_H
#define HPP_GUARD_CTRLPP_PROPAGATE_H

#include "ctrlpp/state_space.h"

namespace ctrlpp {

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
constexpr auto propagate(
    const discrete_state_space<Scalar, NX, NU, NY>& sys,
    const Vector<Scalar, NX>& x,
    const Vector<Scalar, NU>& u)
    -> Vector<Scalar, NX>
{
    return (sys.A * x + sys.B * u).eval();
}

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
constexpr auto output(
    const discrete_state_space<Scalar, NX, NU, NY>& sys,
    const Vector<Scalar, NX>& x,
    const Vector<Scalar, NU>& u)
    -> Vector<Scalar, NY>
{
    return (sys.C * x + sys.D * u).eval();
}

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
constexpr auto output(
    const continuous_state_space<Scalar, NX, NU, NY>& sys,
    const Vector<Scalar, NX>& x,
    const Vector<Scalar, NU>& u)
    -> Vector<Scalar, NY>
{
    return (sys.C * x + sys.D * u).eval();
}

}

#endif
