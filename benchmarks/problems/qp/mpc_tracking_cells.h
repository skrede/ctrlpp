#ifndef HPP_GUARD_BENCHMARKS_PROBLEMS_QP_MPC_TRACKING_CELLS_H
#define HPP_GUARD_BENCHMARKS_PROBLEMS_QP_MPC_TRACKING_CELLS_H

#include "qp/sparse_mpc_sweep.h"

#include <format>
#include <string>
#include <vector>
#include <cstdint>
#include <optional>
#include <string_view>

namespace ctrlpp::bench::problems::qp
{

struct mpc_variant
{
    std::string_view name;
    double input_limit;
    std::optional<double> state_limit;
    double tracking_target;
    double x0_amplitude;
};

inline constexpr double stage_state_weight = 1.0;
inline constexpr double stage_input_weight = 0.1;
inline constexpr double legacy_terminal_weight = 10.0;
inline constexpr std::int32_t long_horizon_min_steps = 50;

inline auto nominal_variant() -> mpc_variant
{
    return {"nominal", 0.5, std::nullopt, 0.0, 5.0};
}

/// The target lies outside the state box, so the optimum holds the actuated states on the box's upper
/// face; every row of the plant has infinity-norm one, so a zero input keeps any initial state inside
/// the box feasible.
inline auto tight_variant() -> mpc_variant
{
    return {"tight", 0.2, 1.0, 2.0, 1.0};
}

inline auto mpc_variants() -> std::vector<mpc_variant>
{
    return {nominal_variant(), tight_variant()};
}

inline auto tracking_cells() -> std::vector<mpc_sweep_cell>
{
    std::vector<mpc_sweep_cell> cells = mpc_sweep_cells();
    cells.insert(cells.end(), {{4, 2, 50}, {4, 2, 100}, {8, 3, 50}, {8, 3, 100}});
    return cells;
}

inline bool is_long_horizon(const mpc_sweep_cell& cell)
{
    return cell.horizon >= long_horizon_min_steps;
}

inline auto cell_id(const mpc_sweep_cell& cell) -> std::string
{
    return std::format("nx{}_nu{}_N{}", cell.nx, cell.nu, cell.horizon);
}

inline auto family_id(const mpc_sweep_cell& cell) -> std::string
{
    return std::format("nx{}_nu{}", cell.nx, cell.nu);
}

}

#endif
