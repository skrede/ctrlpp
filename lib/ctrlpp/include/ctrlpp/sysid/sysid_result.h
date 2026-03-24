#ifndef HPP_GUARD_CTRLPP_SYSID_SYSID_RESULT_H
#define HPP_GUARD_CTRLPP_SYSID_SYSID_RESULT_H

#include "ctrlpp/model/state_space.h"

#include "ctrlpp/sysid/fit_metrics.h"

#include <Eigen/Dense>

namespace ctrlpp {

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct arx_result
{
    discrete_state_space<Scalar, NX, NU, NY> system{};
    fit_metrics<Scalar> metrics{};
};

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct n4sid_result
{
    discrete_state_space<Scalar, NX, NU, NY> system{};
    Eigen::VectorX<Scalar> singular_values{};
    fit_metrics<Scalar> metrics{};
    Scalar condition_number{};
};

}

#endif
