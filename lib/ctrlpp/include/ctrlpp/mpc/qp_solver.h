#ifndef HPP_GUARD_CTRLPP_MPC_QP_SOLVER_H
#define HPP_GUARD_CTRLPP_MPC_QP_SOLVER_H

#include "ctrlpp/mpc/qp_types.h"

#include <concepts>

namespace ctrlpp {

template<typename S>
concept qp_solver = requires {
    typename S::scalar_type;
} && requires(S solver, const qp_problem<typename S::scalar_type>& prob,
              const qp_update<typename S::scalar_type>& upd) {
    { solver.setup(prob) } -> std::same_as<void>;
    { solver.solve(upd) } -> std::same_as<qp_result<typename S::scalar_type>>;
};

}

#endif
