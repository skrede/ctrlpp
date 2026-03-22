#ifndef HPP_GUARD_CTRLPP_MPC_NMPC_CONFIG_H
#define HPP_GUARD_CTRLPP_MPC_NMPC_CONFIG_H

#include "ctrlpp/types.h"

#include <cstddef>
#include <functional>
#include <optional>

namespace ctrlpp {

namespace detail {

template<typename Scalar, std::size_t N>
constexpr auto default_penalty() -> Vector<Scalar, N>
{
    if constexpr (N == 0) {
        return Vector<Scalar, 0>{};
    } else {
        return Vector<Scalar, N>::Constant(Scalar{1e4});
    }
}

}

/// Configuration for nonlinear MPC, mirroring mpc_config with additional
/// custom cost override fields for nonlinear objectives.
template<typename Scalar, std::size_t NX, std::size_t NU,
         std::size_t NC = 0, std::size_t NTC = 0>
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

    /// Path constraint callable g(x, u) -> Vector<NC>, only meaningful when NC > 0.
    /// Constraint is satisfied when g(x, u) <= 0 element-wise.
    std::optional<std::function<Vector<Scalar, NC>(const Vector<Scalar, NX>&,
                                                    const Vector<Scalar, NU>&)>> path_constraint{};

    /// Terminal constraint callable h(x) -> Vector<NTC>, only meaningful when NTC > 0.
    /// Constraint is satisfied when h(x_N) <= 0 element-wise.
    std::optional<std::function<Vector<Scalar, NTC>(const Vector<Scalar, NX>&)>> terminal_constraint{};

    /// When true (default), constraints are softened with slack variables
    /// and L1 penalty to prevent solver failure on infeasible configurations.
    bool soft_constraints{true};

    /// Per-constraint L1 penalty weights for path constraints.
    Vector<Scalar, NC> path_penalty{detail::default_penalty<Scalar, NC>()};

    /// Per-constraint L1 penalty weights for terminal constraints.
    Vector<Scalar, NTC> terminal_penalty{detail::default_penalty<Scalar, NTC>()};
};

}

#endif
