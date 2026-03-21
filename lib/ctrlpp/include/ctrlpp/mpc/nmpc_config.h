#ifndef HPP_GUARD_CTRLPP_MPC_NMPC_CONFIG_H
#define HPP_GUARD_CTRLPP_MPC_NMPC_CONFIG_H

#include "ctrlpp/types.h"

#include <cstddef>
#include <functional>
#include <optional>

namespace ctrlpp {

/// Configuration for nonlinear MPC, mirroring mpc_config with additional
/// custom cost override fields for nonlinear objectives.
template<typename Scalar, std::size_t NX, std::size_t NU>
struct nmpc_config {
    int horizon{1};
    Matrix<Scalar, NX, NX> Q{Matrix<Scalar, NX, NX>::Identity()};
    Matrix<Scalar, NU, NU> R{Matrix<Scalar, NU, NU>::Identity()};
    std::optional<Matrix<Scalar, NX, NX>> Qf{};
    std::optional<Vector<Scalar, NU>> u_min{};
    std::optional<Vector<Scalar, NU>> u_max{};
    std::optional<Vector<Scalar, NX>> x_min{};
    std::optional<Vector<Scalar, NX>> x_max{};
    std::optional<Vector<Scalar, NU>> du_max{};

    /// Custom stage cost override. When set, replaces the quadratic
    /// cost x'Qx + u'Ru entirely for each stage.
    std::optional<std::function<Scalar(const Vector<Scalar, NX>&,
                                       const Vector<Scalar, NU>&)>> stage_cost{};

    /// Custom terminal cost override. When set, replaces the quadratic
    /// terminal cost x'Qf*x entirely.
    std::optional<std::function<Scalar(const Vector<Scalar, NX>&)>> terminal_cost{};
};

}

#endif
