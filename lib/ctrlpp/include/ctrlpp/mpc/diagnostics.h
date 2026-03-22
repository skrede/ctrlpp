#ifndef HPP_GUARD_CTRLPP_MPC_DIAGNOSTICS_H
#define HPP_GUARD_CTRLPP_MPC_DIAGNOSTICS_H

#include "ctrlpp/mpc/qp_types.h"

#include <cstdint>

namespace ctrlpp {

template<typename Scalar>
struct mpc_diagnostics {
    solve_status status{solve_status::error};
    int iterations{};
    Scalar solve_time{};
    Scalar cost;
    Scalar primal_residual;
    Scalar dual_residual;
    Scalar max_constraint_violation;
    Scalar max_path_constraint_violation{};
    Scalar max_terminal_constraint_violation{};
    Scalar total_slack{};
};

}

#endif
