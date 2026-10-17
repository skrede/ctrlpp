#ifndef HPP_GUARD_CTRLPP_MPC_OSQP_SOLVER_H
#define HPP_GUARD_CTRLPP_MPC_OSQP_SOLVER_H

/// @brief OSQP-backed QP solver wrapper for linear MPC.
///
/// @cite stellato2020 -- Stellato, Banjac, Goulart, Bemporad & Boyd, "OSQP: An Operator Splitting Solver for Quadratic Programs", Math. Prog. Comp. 12(4):637-672, 2020

#include "ctrlpp/expected.h"

#include "ctrlpp/mpc/qp_types.h"

#include "ctrlpp/mpc/detail/osqp_interop.h"

#include <Eigen/Sparse>

#include <osqp.h>

#include <cstdint>

namespace ctrlpp
{

class osqp_solver
{
public:
    using scalar_type = double;

    explicit osqp_solver(double eps_abs = 1e-3, double eps_rel = 1e-3, std::int32_t max_iter = 4000, bool verbose = false, bool warm_starting = true, bool polishing = true)
        : p_upper_{}
        , solver_{}
        , a_{}
        , max_iter_{max_iter}
        , warm_starting_{warm_starting}
        , eps_abs_{eps_abs}
        , eps_rel_{eps_rel}
        , polishing_{polishing}
        , verbose_{verbose}
        , n_{0}
        , m_{0}
    {
    }

    /// @brief Preset constructor: `accuracy` polishes, `speed` does not. See
    /// `qp_preset`. Backend-agnostic with `argmin_qp_solver`, so a call site can
    /// pick the accuracy/speed tradeoff independently of the chosen backend.
    explicit osqp_solver(qp_preset preset)
        : osqp_solver(1e-3, 1e-3, 4000, false, true, preset == qp_preset::accuracy)
    {
    }

    osqp_solver(const osqp_solver&) = delete;
    auto operator=(const osqp_solver&) -> osqp_solver& = delete;

    osqp_solver(osqp_solver&&) noexcept = default;
    auto operator=(osqp_solver&&) noexcept -> osqp_solver& = default;

    /// @brief Fallible setup: initializes the OSQP workspace from the problem data.
    /// Returns an empty expected on success and an `osqp_setup_error` on failure.
    /// This is the only setup shape, and it is available in every build mode,
    /// including `-fno-exceptions`.
    auto setup(const qp_problem<double>& problem) -> ctrlpp::expected<void, osqp_setup_error>
    {
        solver_.reset();
        prepare_sparse_matrices(problem);
        return configure_and_create_solver(problem);
    }

    auto solve(const qp_update<double>& update) -> qp_result<double>
    {
        osqp_update_data_vec(solver_.get(), update.q.size() > 0 ? update.q.data() : nullptr, update.l.size() > 0 ? update.l.data() : nullptr, update.u.size() > 0 ? update.u.data() : nullptr);

        if(update.warm_x.size() > 0 || update.warm_y.size() > 0)
        {
            osqp_warm_start(solver_.get(), update.warm_x.size() > 0 ? update.warm_x.data() : nullptr, update.warm_y.size() > 0 ? update.warm_y.data() : nullptr);
        }

        osqp_solve(solver_.get());
        return detail::read_osqp_result(*solver_, n_, m_);
    }

private:
    detail::osqp_sparse_matrix p_upper_;
    detail::osqp_workspace solver_;
    detail::osqp_sparse_matrix a_;
    std::int32_t max_iter_;
    bool warm_starting_;
    double eps_abs_;
    double eps_rel_;
    bool polishing_;
    bool verbose_;
    OSQPInt n_;
    OSQPInt m_;

    void prepare_sparse_matrices(const qp_problem<double>& problem)
    {
        p_upper_ = problem.P.template triangularView<Eigen::Upper>();
        p_upper_.makeCompressed();

        a_ = problem.A.template cast<OSQPFloat>();
        a_.makeCompressed();

        n_ = static_cast<OSQPInt>(problem.P.cols());
        m_ = static_cast<OSQPInt>(problem.A.rows());
    }

    auto make_settings() const -> detail::osqp_settings
    {
        detail::osqp_settings settings{OSQPSettings_new()};
        if(settings)
        {
            settings->eps_abs = eps_abs_;
            settings->eps_rel = eps_rel_;
            settings->max_iter = max_iter_;
            settings->verbose = verbose_;
            settings->warm_starting = warm_starting_;
            settings->polishing = polishing_;
        }
        return settings;
    }

    auto configure_and_create_solver(const qp_problem<double>& problem) -> ctrlpp::expected<void, osqp_setup_error>
    {
        auto p_csc = detail::make_csc(p_upper_);
        auto a_csc = detail::make_csc(a_);

        const auto settings = make_settings();
        if(!settings)
            return ctrlpp::unexpected(osqp_setup_error::settings_allocation_failed);

        OSQPSolver* solver = nullptr;
        auto exit_flag = osqp_setup(&solver, &p_csc, const_cast<OSQPFloat*>(problem.q.data()), &a_csc, const_cast<OSQPFloat*>(problem.l.data()), const_cast<OSQPFloat*>(problem.u.data()), m_, n_, settings.get());
        solver_.reset(solver);

        if(exit_flag != 0)
        {
            // osqp_setup publishes the solver pointer before the allocations that can fail, and its error paths free nothing,
            // so a failed setup hands back a partially built solver that is ours to release. Dropping it instead would leak it.
            solver_.reset();
            return ctrlpp::unexpected(osqp_setup_error::setup_failed);
        }

        return {};
    }
};

}

#endif
