#ifndef HPP_GUARD_CTRLPP_MHE_MHE_DIAGNOSTICS_H
#define HPP_GUARD_CTRLPP_MHE_MHE_DIAGNOSTICS_H

#include "ctrlpp/mpc/qp_types.h"

#include <cstdint>

namespace ctrlpp
{

template <typename Scalar>
struct mhe_diagnostics
{
    mhe_diagnostics()
        : status{solve_status::error}
        , iterations{}
        , solve_time{}
        , cost{}
        , primal_residual{}
        , dual_residual{}
        , max_constraint_violation{}
        , max_residual_bound_violation{}
        , total_slack{}
        , used_ekf_fallback{false}
    {
    }

    solve_status status;
    std::int32_t iterations;
    Scalar solve_time;
    Scalar cost;
    Scalar primal_residual;
    Scalar dual_residual;
    Scalar max_constraint_violation;
    Scalar max_residual_bound_violation;
    Scalar total_slack;
    bool used_ekf_fallback;
};

namespace detail
{

template <typename Scalar>
mhe_diagnostics<Scalar> ekf_fallback_diagnostics(solve_status status)
{
    mhe_diagnostics<Scalar> diagnostics;
    diagnostics.status = status;
    diagnostics.used_ekf_fallback = true;
    return diagnostics;
}

}

}

#endif
