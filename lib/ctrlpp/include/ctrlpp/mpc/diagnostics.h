#ifndef HPP_GUARD_CTRLPP_MPC_DIAGNOSTICS_H
#define HPP_GUARD_CTRLPP_MPC_DIAGNOSTICS_H

#include "ctrlpp/mpc/qp_types.h"

#include <cstdint>

namespace ctrlpp
{

template <typename Scalar>
struct mpc_diagnostics
{
    solve_status status{solve_status::error};
    int iterations{};
    Scalar solve_time{};
    Scalar cost{};
    Scalar primal_residual{};
    Scalar dual_residual{};
    Scalar max_constraint_violation{};
    Scalar max_path_constraint_violation{};
    Scalar max_terminal_constraint_violation{};
    Scalar total_slack{};
    /// Set when no terminal weight was configured AND the Riccati solve that
    /// would have supplied one did not produce a solution, so the state weight
    /// stands in for the terminal cost. Every solve after that optimizes a
    /// different problem than an infinite-horizon terminal cost poses, and any
    /// stability argument resting on that terminal cost no longer holds.
    /// Latched at construction, so it reads the same on every solve.
    bool used_state_weight_terminal_cost{false};
};

}

#endif
