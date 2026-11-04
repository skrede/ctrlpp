#ifndef HPP_GUARD_CTRLPP_MPC_MPC_CONFIG_H
#define HPP_GUARD_CTRLPP_MPC_MPC_CONFIG_H

#include "ctrlpp/mpc/terminal_set.h"

#include "ctrlpp/types.h"

#include <cstddef>
#include <cstdint>
#include <optional>

namespace ctrlpp
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY = NX>
struct mpc_config
{
    std::int32_t horizon{1};
    Matrix<Scalar, NY, NY> Q{Matrix<Scalar, NY, NY>::Identity()};
    Matrix<Scalar, NU, NU> R{Matrix<Scalar, NU, NU>::Identity()};
    std::optional<Matrix<Scalar, NY, NY>> Qf{};
    std::optional<Vector<Scalar, NU>> u_min{};
    std::optional<Vector<Scalar, NU>> u_max{};
    std::optional<Vector<Scalar, NX>> x_min{};
    std::optional<Vector<Scalar, NX>> x_max{};
    std::optional<Vector<Scalar, NU>> du_max{};
    Scalar soft_penalty{Scalar{1e4}};
    std::optional<Vector<Scalar, NX>> soft_state_penalty{};
    std::optional<terminal_set<Scalar, NX>> terminal_constraint_set{};
    bool hard_state_constraints{false};
};

}

#endif
