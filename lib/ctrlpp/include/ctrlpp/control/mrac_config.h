#ifndef HPP_GUARD_CTRLPP_CONTROL_MRAC_CONFIG_H
#define HPP_GUARD_CTRLPP_CONTROL_MRAC_CONFIG_H

#include "ctrlpp/types.h"

#include "ctrlpp/control/mrac_policies.h"

#include "ctrlpp/model/state_space.h"

#include <cstddef>
#include <type_traits>

namespace ctrlpp
{

namespace detail
{

template <typename Robustification, typename Scalar, typename = void>
struct robustification_options;

template <typename R, typename Scalar>
struct robustification_options<R, Scalar, std::void_t<typename R::options_t>>
{
    using type = typename R::options_t;
};

template <typename R, typename Scalar>
struct robustification_options<R, Scalar, std::void_t<typename R::template options_t<Scalar>>>
{
    using type = typename R::template options_t<Scalar>;
};

template <typename R, typename Scalar>
using robustification_options_t = typename robustification_options<R, Scalar>::type;

}

template <typename Scalar, std::size_t NX, std::size_t NU,
          typename Robustification = no_robustification>
struct mrac_config
{
    discrete_state_space<Scalar, NX, NU, NX> reference_model{};
    Scalar gamma{};
    Scalar sign_b{1};
    Scalar theta_x_0{};
    Scalar theta_r_0{};
    Vector<Scalar, NX> x_model_0 = Vector<Scalar, NX>::Zero();
    detail::robustification_options_t<Robustification, Scalar> robustification{};
};

}

#endif
