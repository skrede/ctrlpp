#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_OSQP_INTEROP_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_OSQP_INTEROP_H

#include "ctrlpp/mpc/qp_types.h"

#include <Eigen/Sparse>

#include <osqp.h>

#include <memory>
#include <type_traits>

namespace ctrlpp::detail
{

static_assert(std::is_same_v<double, OSQPFloat>, "OSQP must be compiled with double precision");

using osqp_sparse_matrix = Eigen::SparseMatrix<OSQPFloat, Eigen::ColMajor, OSQPInt>;

struct osqp_workspace_deleter
{
    void operator()(OSQPSolver* solver) const noexcept { osqp_cleanup(solver); }
};

struct osqp_settings_deleter
{
    void operator()(OSQPSettings* settings) const noexcept { OSQPSettings_free(settings); }
};

using osqp_workspace = std::unique_ptr<OSQPSolver, osqp_workspace_deleter>;
using osqp_settings = std::unique_ptr<OSQPSettings, osqp_settings_deleter>;

inline auto make_csc(const osqp_sparse_matrix& mat) -> OSQPCscMatrix
{
    OSQPCscMatrix csc{};
    csc.m = static_cast<OSQPInt>(mat.rows());
    csc.n = static_cast<OSQPInt>(mat.cols());
    csc.nzmax = static_cast<OSQPInt>(mat.nonZeros());
    csc.x = const_cast<OSQPFloat*>(mat.valuePtr());
    csc.i = const_cast<OSQPInt*>(mat.innerIndexPtr());
    csc.p = const_cast<OSQPInt*>(mat.outerIndexPtr());
    csc.nz = -1;
    return csc;
}

constexpr auto map_osqp_status(OSQPInt status_val) -> solve_status
{
    switch(status_val)
    {
    case OSQP_SOLVED:
        return solve_status::optimal;
    case OSQP_SOLVED_INACCURATE:
        return solve_status::solved_inaccurate;
    case OSQP_PRIMAL_INFEASIBLE:
        [[fallthrough]];
    case OSQP_PRIMAL_INFEASIBLE_INACCURATE:
        return solve_status::infeasible;
    case OSQP_DUAL_INFEASIBLE:
        [[fallthrough]];
    case OSQP_DUAL_INFEASIBLE_INACCURATE:
        return solve_status::unbounded;
    case OSQP_MAX_ITER_REACHED:
        return solve_status::max_iterations;
    case OSQP_TIME_LIMIT_REACHED:
        return solve_status::time_limit;
    case OSQP_NON_CVX:
        return solve_status::non_convex;
    default:
        return solve_status::error;
    }
}

inline auto read_osqp_result(const OSQPSolver& solver, OSQPInt n, OSQPInt m) -> qp_result<double>
{
    qp_result<double> result;
    result.x = Eigen::Map<const Eigen::VectorXd>(solver.solution->x, n);
    result.y = Eigen::Map<const Eigen::VectorXd>(solver.solution->y, m);
    result.status = map_osqp_status(solver.info->status_val);
    result.objective = solver.info->obj_val;
    result.solve_time = solver.info->solve_time;
    result.iterations = static_cast<int>(solver.info->iter);
    result.primal_residual = solver.info->prim_res;
    result.dual_residual = solver.info->dual_res;
    return result;
}

}

#endif
