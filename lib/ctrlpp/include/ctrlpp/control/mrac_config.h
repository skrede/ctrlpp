#ifndef HPP_GUARD_CTRLPP_CONTROL_MRAC_CONFIG_H
#define HPP_GUARD_CTRLPP_CONTROL_MRAC_CONFIG_H

#include "ctrlpp/types.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/control/mrac_policies.h"

#include "ctrlpp/model/state_space.h"

#include <cstddef>
#include <type_traits>

namespace ctrlpp
{

namespace detail
{

template <typename R, typename Scalar>
concept has_scalar_options = requires { typename R::template options_t<Scalar>; };

template <typename R, typename Scalar>
consteval auto deduce_robustification_options()
{
    if constexpr (has_scalar_options<R, Scalar>)
        return std::type_identity<typename R::template options_t<Scalar>>{};
    else
        return std::type_identity<typename R::options_t>{};
}

template <typename R, typename Scalar>
using robustification_options_t = typename decltype(deduce_robustification_options<R, Scalar>())::type;

}

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU,
          typename Robustification = no_robustification>
struct mrac_config
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    discrete_state_space<Scalar, NX, NU, NX> reference_model{};
    Matrix<Scalar, NX, NX> gamma_x = Matrix<Scalar, NX, NX>::Zero();
    Matrix<Scalar, NU, NU> gamma_r = Matrix<Scalar, NU, NU>::Zero();
    Matrix<Scalar, NU, NU> sign_b = Matrix<Scalar, NU, NU>::Identity();
    Matrix<Scalar, NU, NX> theta_x_0 = Matrix<Scalar, NU, NX>::Zero();
    Matrix<Scalar, NU, NU> theta_r_0 = Matrix<Scalar, NU, NU>::Zero();
    Vector<Scalar, NX> x_model_0 = Vector<Scalar, NX>::Zero();
    Matrix<Scalar, NX, NX> W = Matrix<Scalar, NX, NX>::Identity();
    detail::robustification_options_t<Robustification, Scalar> robustification{};
};

}

#endif
