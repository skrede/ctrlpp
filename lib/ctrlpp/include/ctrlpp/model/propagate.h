#ifndef HPP_GUARD_CTRLPP_MODEL_PROPAGATE_H
#define HPP_GUARD_CTRLPP_MODEL_PROPAGATE_H

/// @brief Discrete state-space propagation and output-equation evaluation.
///
/// @cite kailath1980 -- Kailath, "Linear Systems", 1980, Sec. 2.1 (linear state-space evolution)
/// @cite franklin2015 -- Franklin et al., "Feedback Control of Dynamic Systems", 2015

#include "ctrlpp/model/state_space.h"

namespace ctrlpp
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
constexpr Vector<Scalar, NX> propagate(const discrete_state_space<Scalar, NX, NU, NY>& sys, const Vector<Scalar, NX>& x, const Vector<Scalar, NU>& u)
{
    return (sys.A * x + sys.B * u).eval();
}

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
constexpr Vector<Scalar, NY> output(const discrete_state_space<Scalar, NX, NU, NY>& sys, const Vector<Scalar, NX>& x, const Vector<Scalar, NU>& u)
{
    return (sys.C * x + sys.D * u).eval();
}

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
constexpr Vector<Scalar, NY> output(const continuous_state_space<Scalar, NX, NU, NY>& sys, const Vector<Scalar, NX>& x, const Vector<Scalar, NU>& u)
{
    return (sys.C * x + sys.D * u).eval();
}

}

#endif
