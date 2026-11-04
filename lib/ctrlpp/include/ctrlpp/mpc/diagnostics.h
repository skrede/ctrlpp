#ifndef HPP_GUARD_CTRLPP_MPC_DIAGNOSTICS_H
#define HPP_GUARD_CTRLPP_MPC_DIAGNOSTICS_H

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/nlp_types.h"

#include <cstdint>

namespace ctrlpp
{

template <typename Scalar>
struct mpc_diagnostics
{
    mpc_diagnostics()
        : status{solve_status::error}
        , iterations{}
        , solve_time{}
        , cost{}
        , primal_residual{}
        , dual_residual{}
        , max_constraint_violation{}
        , max_path_constraint_violation{}
        , max_terminal_constraint_violation{}
        , total_slack{}
        , stop_criterion{nlp_stop_criterion::unreported}
    {
    }

    solve_status status;
    std::int32_t iterations;
    Scalar solve_time;
    Scalar cost;
    Scalar primal_residual;
    Scalar dual_residual;
    Scalar max_constraint_violation;
    Scalar max_path_constraint_violation;
    Scalar max_terminal_constraint_violation;
    Scalar total_slack;
    // Stays `unreported` for the QP-backed controllers and for nonlinear backends that do not say which test stopped.
    nlp_stop_criterion stop_criterion;
};

}

#endif
