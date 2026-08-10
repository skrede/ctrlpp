#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_ARGMIN_SHOOTING_JACOBIAN_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_ARGMIN_SHOOTING_JACOBIAN_H

#include "ctrlpp/types.h"

#include <Eigen/Core>

#include <span>
#include <cstddef>
#include <algorithm>
#include <functional>

namespace ctrlpp::bench
{

/// Where each block of the multiple-shooting decision vector and of its
/// continuity constraints starts. The states occupy the leading block, so their
/// offset is zero and is not carried.
struct shooting_layout
{
    int nu;
    int nx;
    int n_vars;
    int horizon;
    int u_offset;
    int n_constraints;
};

template <std::size_t NX, std::size_t NU>
auto shooting_layout_of(int horizon) -> shooting_layout
{
    const int nx = static_cast<int>(NX);
    const int nu = static_cast<int>(NU);
    return {nu, nx, (horizon + 1) * nx + horizon * nu, horizon, (horizon + 1) * nx, (horizon + 1) * nx};
}

/// The flattened Jacobian is column-major over (n_constraints x n_vars), which
/// is the layout the solver bridge reads it back in.
inline auto entry_of(const shooting_layout& layout, int row, int column) -> std::size_t
{
    return static_cast<std::size_t>(row) + static_cast<std::size_t>(column)
                                               * static_cast<std::size_t>(layout.n_constraints);
}

inline void write_initial_state_block(const shooting_layout& layout, std::span<double> jacobian)
{
    for(int i = 0; i < layout.nx; ++i)
        jacobian[entry_of(layout, i, i)] = 1.0;
}

template <std::size_t NX, std::size_t NU, typename Dynamics>
void write_continuity_block(const shooting_layout& layout, const Dynamics& dynamics, int k,
                            std::span<const double> z, std::span<double> jacobian)
{
    const Eigen::Map<const ctrlpp::Vector<double, NX>> xk(z.data() + k * layout.nx);
    const Eigen::Map<const ctrlpp::Vector<double, NU>> uk(z.data() + layout.u_offset + k * layout.nu);
    const auto state_partials = dynamics.jacobian_x(xk, uk);
    const auto input_partials = dynamics.jacobian_u(xk, uk);

    for(int i = 0; i < layout.nx; ++i)
    {
        const int row = (k + 1) * layout.nx + i;
        for(int j = 0; j < layout.nx; ++j)
            jacobian[entry_of(layout, row, k * layout.nx + j)] = -state_partials(i, j);
        jacobian[entry_of(layout, row, (k + 1) * layout.nx + i)] = 1.0;
        for(int j = 0; j < layout.nu; ++j)
            jacobian[entry_of(layout, row, layout.u_offset + k * layout.nu + j)] = -input_partials(i, j);
    }
}

/// The exact partials of the initial-state and continuity constraints, handed to
/// the solver in place of the finite differences it would otherwise form.
template <std::size_t NX, std::size_t NU, typename Dynamics>
auto build_shooting_jacobian(const Dynamics& dynamics, int horizon)
    -> std::function<void(std::span<const double>, std::span<double>)>
{
    const shooting_layout layout = shooting_layout_of<NX, NU>(horizon);
    return [=](std::span<const double> z, std::span<double> jacobian)
    {
        std::fill(jacobian.begin(), jacobian.end(), 0.0);
        write_initial_state_block(layout, jacobian);
        for(int k = 0; k < layout.horizon; ++k)
            write_continuity_block<NX, NU>(layout, dynamics, k, z, jacobian);
    };
}

}

#endif
