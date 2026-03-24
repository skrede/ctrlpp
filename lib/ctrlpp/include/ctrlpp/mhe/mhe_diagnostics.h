#ifndef HPP_GUARD_CTRLPP_MHE_MHE_DIAGNOSTICS_H
#define HPP_GUARD_CTRLPP_MHE_MHE_DIAGNOSTICS_H

#include "ctrlpp/mpc/qp_types.h"

namespace ctrlpp
{

/// Diagnostics for MHE solves, mirroring mpc_diagnostics with MHE-specific fields.
template <typename Scalar>
struct mhe_diagnostics
{
    solve_status status{solve_status::error};
    int iterations{};
    Scalar solve_time{};
    Scalar cost{};
    Scalar primal_residual{};
    Scalar dual_residual{};
    Scalar max_constraint_violation{};
    Scalar max_residual_bound_violation{};
    Scalar total_slack{};
    bool used_ekf_fallback{false};
};

} // namespace ctrlpp

#endif
