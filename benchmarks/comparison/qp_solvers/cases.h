#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_CASES_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_CASES_H

#include "instance.h"
#include "registry.h"

#include "qp/mpc_tracking_cells.h"

#include <Eigen/Dense>

#include <memory>
#include <string>
#include <cstdint>

namespace ctrlpp::bench::qp_solvers
{

/// A problem's inequality rows in the uncondensed form shared by both cohorts: the input boxes, then
/// the state boxes, with the state rows' values and bounds in state coordinates. Binding is judged on
/// these, so an instance binds identically whichever cohort poses it.
struct inequality_rows
{
    Eigen::VectorXd value;
    Eigen::VectorXd lower;
    Eigen::VectorXd upper;
};

class problem_case
{
public:
    virtual ~problem_case() = default;

    virtual auto instance() const -> const qp_instance& = 0;

    virtual void fill(const Eigen::VectorXd& x0, qp_tick& tick) const = 0;

    virtual auto rows_of(const qp_tick& tick, const Eigen::VectorXd& x) const -> inequality_rows = 0;

    virtual auto state_dimension() const -> std::int32_t = 0;

    virtual auto amplitude() const -> double = 0;
};

auto make_sparse_case(const problems::qp::mpc_sweep_cell& cell, const problems::qp::mpc_variant& variant)
    -> std::unique_ptr<problem_case>;

auto make_condensed_case(const problems::qp::mpc_sweep_cell& cell, const problems::qp::mpc_variant& variant)
    -> std::unique_ptr<problem_case>;

/// dense_mpc_shaped: no initial state, every tick identical.
auto make_shaped_case() -> std::unique_ptr<problem_case>;

auto initial_state(const problem_case& source, regime kind, std::int32_t seed, std::int32_t step) -> Eigen::VectorXd;

inline constexpr std::string_view shaped_cell_id = "dense_mpc_shaped";

}

#endif
