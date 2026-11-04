#ifndef HPP_GUARD_CTRLPP_MPC_ARGMIN_QP_SOLVER_H
#define HPP_GUARD_CTRLPP_MPC_ARGMIN_QP_SOLVER_H

#include "ctrlpp/config.h"

#if defined(CTRLPP_HAS_ARGMIN) && __has_include(<argmin/qp/sparse_admm_qp.h>)

    #define CTRLPP_HAS_ARGMIN_QP 1

    #include "ctrlpp/mpc/qp_types.h"

    #include "ctrlpp/expected.h"

    #include <argmin/qp/qp_types.h>
    #include <argmin/qp/sparse_admm_qp.h>

    #include <argmin/options/sparse_qp_options.h>

    #include <cstdint>

namespace ctrlpp
{

enum class argmin_qp_setup_error : std::uint8_t
{
    pose_failed,
    invalid_options,
};

// @cite stellato2020 -- Stellato et al., "OSQP: An Operator Splitting Solver for Quadratic Programs", Math. Prog. Comp. 12(4), 2020
class argmin_qp_solver
{
public:
    using scalar_type = double;

    // This binds only argmin's stable QP contract (tolerances, budget, warm start, polish); the operator-splitting
    // knobs are the volatile surface argmin intends to demote, so a reshape of them is a no-op here.
    // Every negative std::int32_t converts to a value above INT32_MAX, which argmin's option-domain check refuses,
    // so the conversion never wraps a budget into the valid range.
    explicit argmin_qp_solver(double eps_abs = 1e-3, double eps_rel = 1e-3, std::int32_t max_iter = 4000, bool /*verbose*/ = false, bool warm_starting = true, bool polishing = true)
        : m_solver{}
        , m_opts{.eps_abs = eps_abs, .eps_rel = eps_rel, .max_iterations = static_cast<std::uint32_t>(max_iter), .polish = polishing, .warm_start = warm_starting}
        , m_last{}
        , m_posed{false}
    {
    }

    explicit argmin_qp_solver(qp_preset preset)
        : argmin_qp_solver(1e-3, 1e-3, 4000, false, true, preset == qp_preset::accuracy)
    {
    }

    ctrlpp::expected<void, argmin_qp_setup_error> setup(const qp_problem<double>& problem)
    {
        // qp_problem carries P upper-triangular, OSQP's convention; argmin reads P as full symmetric.
        const Eigen::SparseMatrix<double> full_p = problem.P.template selfadjointView<Eigen::Upper>();
        if(auto err = m_solver.solve_into(full_p, problem.q, problem.A, problem.l, problem.u, m_last, m_opts))
            return ctrlpp::unexpected(translate_error(*err));
        m_posed = true;
        return {};
    }

    qp_result<double> solve(const qp_update<double>& update)
    {
        if(!m_posed)
            return {};

        if(update.warm_x.size() > 0 && update.warm_y.size() > 0)
            m_solver.warm_start(update.warm_x, update.warm_y);

        if(auto err = m_solver.resolve_into(update.q, update.l, update.u, m_last, m_opts))
            return {};

        return translate(m_last);
    }

private:
    argmin::sparse_admm_qp_solver<double> m_solver;
    argmin::sparse_qp_options<> m_opts;
    argmin::qp_result<double> m_last;
    bool m_posed;

    static qp_result<double> translate(const argmin::qp_result<double>& r)
    {
        qp_result<double> out;
        out.status = translate_status(r.status);
        out.x = r.x;
        out.y = r.y;
        out.objective = r.objective_value;
        // argmin takes no timing measurements.
        out.solve_time = double{0};
        out.iterations = r.iterations;
        out.primal_residual = r.primal_residual;
        out.dual_residual = r.dual_residual;
        return out;
    }

    static solve_status translate_status(argmin::qp_solve_status s)
    {
        switch(s)
        {
            case argmin::qp_solve_status::solved:
                return solve_status::optimal;
            case argmin::qp_solve_status::solved_inaccurate:
                return solve_status::solved_inaccurate;
            case argmin::qp_solve_status::max_iterations:
                return solve_status::max_iterations;
            case argmin::qp_solve_status::primal_infeasible:
                return solve_status::infeasible;
            case argmin::qp_solve_status::dual_infeasible:
                return solve_status::unbounded;
        }
        return solve_status::error;
    }

    static argmin_qp_setup_error translate_error(argmin::qp_error e)
    {
        switch(e)
        {
            case argmin::qp_error::invalid_options:
                return argmin_qp_setup_error::invalid_options;
            case argmin::qp_error::dimension_mismatch:
            case argmin::qp_error::invalid_bounds:
            case argmin::qp_error::non_finite_input:
            case argmin::qp_error::capacity_exceeded:
            case argmin::qp_error::infeasible_start:
            case argmin::qp_error::invalid_problem:
            case argmin::qp_error::factorization_failed:
                return argmin_qp_setup_error::pose_failed;
        }
        return argmin_qp_setup_error::pose_failed;
    }
};

}

#endif // CTRLPP_HAS_ARGMIN && __has_include(<argmin/qp/sparse_admm_qp.h>)

#endif
