#ifndef HPP_GUARD_CTRLPP_SYSID_SYSID_RESULT_H
#define HPP_GUARD_CTRLPP_SYSID_SYSID_RESULT_H

#include "ctrlpp/sysid/fit_metrics.h"

#include "ctrlpp/types.h"

#include "ctrlpp/model/state_space.h"

#include <Eigen/Dense>

namespace ctrlpp
{

template <typename Scalar>
struct arx_diagnostics
{
    arx_diagnostics()
        : effective_samples{}
        , parameter_count{}
        , numerical_rank{}
        , residual_norm{}
    {
    }

    std::size_t effective_samples;
    std::size_t parameter_count;
    std::size_t numerical_rank;
    Scalar residual_norm;
};

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct arx_result
{
    arx_result()
        : system{}
        , metrics{}
        , diagnostics{}
    {
    }

    discrete_state_space<Scalar, NX, NU, NY> system;
    fit_metrics<Scalar> metrics;
    arx_diagnostics<Scalar> diagnostics;
};

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct moesp_result
{
    moesp_result()
        : system{}
        , singular_values{}
        , metrics{}
        , condition_number{}
    {
    }

    discrete_state_space<Scalar, NX, NU, NY> system;
    Eigen::VectorX<Scalar> singular_values;
    fit_metrics<Scalar> metrics;
    Scalar condition_number;
    Vector<Scalar, NX> initial_state;
};

}

#endif
