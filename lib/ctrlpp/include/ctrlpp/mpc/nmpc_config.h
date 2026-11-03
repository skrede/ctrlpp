#ifndef HPP_GUARD_CTRLPP_MPC_NMPC_CONFIG_H
#define HPP_GUARD_CTRLPP_MPC_NMPC_CONFIG_H

#include "ctrlpp/types.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <functional>

namespace ctrlpp
{

namespace detail
{

template <typename Scalar, std::size_t N>
inline Vector<Scalar, N> default_penalty()
{
    if constexpr(N == 0)
    {
        return Vector<Scalar, 0>{};
    }
    else
    {
        return Vector<Scalar, N>::Constant(Scalar{1e4});
    }
}

}

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC = 0, std::size_t NTC = 0>
struct nmpc_config
{
    std::int32_t horizon{1};
    Matrix<Scalar, NX, NX> Q{Matrix<Scalar, NX, NX>::Identity()};
    Matrix<Scalar, NU, NU> R{Matrix<Scalar, NU, NU>::Identity()};
    std::optional<Matrix<Scalar, NX, NX>> Qf{};
    std::optional<Vector<Scalar, NU>> u_min{};
    std::optional<Vector<Scalar, NU>> u_max{};
    std::optional<Vector<Scalar, NX>> x_min{};
    std::optional<Vector<Scalar, NX>> x_max{};
    std::optional<Vector<Scalar, NU>> du_max{};

    // A set cost callback replaces the quadratic stage or terminal cost entirely.
    // A constraint callback is satisfied where it is <= 0 element-wise.
    std::optional<std::function<Scalar(const Vector<Scalar, NX>&, const Vector<Scalar, NU>&)>> stage_cost_cb{};
    std::optional<std::function<Scalar(const Vector<Scalar, NX>&)>> terminal_cost_cb{};
    std::optional<std::function<Vector<Scalar, NC>(const Vector<Scalar, NX>&, const Vector<Scalar, NU>&)>> path_constraint_cb{};
    std::optional<std::function<Vector<Scalar, NTC>(const Vector<Scalar, NX>&)>> terminal_constraint_cb{};

    // Softens the path and terminal constraints with L1-penalized slack variables.
    bool soft_constraints{true};
    Vector<Scalar, NC> path_penalty{detail::default_penalty<Scalar, NC>()};
    Vector<Scalar, NTC> terminal_penalty{detail::default_penalty<Scalar, NTC>()};
};

}

#endif
