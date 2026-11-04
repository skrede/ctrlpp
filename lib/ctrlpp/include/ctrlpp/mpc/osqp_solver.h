#ifndef HPP_GUARD_CTRLPP_MPC_OSQP_SOLVER_H
#define HPP_GUARD_CTRLPP_MPC_OSQP_SOLVER_H

// @cite stellato2020 -- Stellato, Banjac, Goulart, Bemporad & Boyd, "OSQP: An Operator Splitting Solver for Quadratic Programs", Math. Prog. Comp. 12(4):637-672, 2020

#include "ctrlpp/mpc/qp_types.h"

#include "ctrlpp/mpc/detail/osqp_interop.h"

#include "ctrlpp/expected.h"

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
        : m_p_upper{}
        , m_solver{}
        , m_a{}
        , m_osqp_lower{}
        , m_osqp_upper{}
        , m_max_iter{max_iter}
        , m_warm_starting{warm_starting}
        , m_eps_abs{eps_abs}
        , m_eps_rel{eps_rel}
        , m_polishing{polishing}
        , m_verbose{verbose}
        , m_n_vars{0}
        , m_n_constraints{0}
    {
    }

    explicit osqp_solver(qp_preset preset)
        : osqp_solver(1e-3, 1e-3, 4000, false, true, preset == qp_preset::accuracy)
    {
    }

    osqp_solver(const osqp_solver&) = delete;
    osqp_solver& operator=(const osqp_solver&) = delete;

    osqp_solver(osqp_solver&&) noexcept = default;
    osqp_solver& operator=(osqp_solver&&) noexcept = default;

    ctrlpp::expected<void, osqp_setup_error> setup(const qp_problem<double>& problem)
    {
        m_solver.reset();
        prepare_sparse_matrices(problem);
        if(problem.l.size() != m_n_constraints || problem.u.size() != m_n_constraints)
            return ctrlpp::unexpected(osqp_setup_error::setup_failed);
        m_osqp_lower.resize(m_n_constraints);
        m_osqp_upper.resize(m_n_constraints);
        detail::to_osqp_lower(problem.l, m_osqp_lower);
        detail::to_osqp_upper(problem.u, m_osqp_upper);
        return configure_and_create_solver(problem);
    }

    qp_result<double> solve(const qp_update<double>& update)
    {
        if(!m_solver || !accepts_bound(update.l) || !accepts_bound(update.u))
            return {};

        // OSQP refuses l > u before touching q, l or u, so a refusal left unchecked would solve the previous problem.
        const OSQPFloat* q = update.q.size() > 0 ? update.q.data() : nullptr;
        if(osqp_update_data_vec(m_solver.get(), q, stage(update.l, m_osqp_lower, detail::to_osqp_lower), stage(update.u, m_osqp_upper, detail::to_osqp_upper)) != 0)
            return {};

        if(update.warm_x.size() > 0 || update.warm_y.size() > 0)
            osqp_warm_start(m_solver.get(), update.warm_x.size() > 0 ? update.warm_x.data() : nullptr, update.warm_y.size() > 0 ? update.warm_y.data() : nullptr);

        osqp_solve(m_solver.get());
        return detail::read_osqp_result(*m_solver, m_n_vars, m_n_constraints);
    }

private:
    detail::osqp_sparse_matrix m_p_upper;
    detail::osqp_workspace m_solver;
    detail::osqp_sparse_matrix m_a;
    Eigen::VectorXd m_osqp_lower;
    Eigen::VectorXd m_osqp_upper;
    std::int32_t m_max_iter;
    bool m_warm_starting;
    double m_eps_abs;
    double m_eps_rel;
    bool m_polishing;
    bool m_verbose;
    OSQPInt m_n_vars;
    OSQPInt m_n_constraints;

    void prepare_sparse_matrices(const qp_problem<double>& problem)
    {
        m_p_upper = problem.P.template triangularView<Eigen::Upper>();
        m_p_upper.makeCompressed();

        m_a = problem.A.template cast<OSQPFloat>();
        m_a.makeCompressed();

        m_n_vars = static_cast<OSQPInt>(problem.P.cols());
        m_n_constraints = static_cast<OSQPInt>(problem.A.rows());
    }

    bool accepts_bound(const Eigen::VectorXd& bound) const { return bound.size() == 0 || bound.size() == m_n_constraints; }

    static const OSQPFloat* stage(const Eigen::VectorXd& bound, Eigen::VectorXd& buffer, void (&map)(const Eigen::VectorXd&, Eigen::VectorXd&))
    {
        if(bound.size() == 0)
            return nullptr;
        map(bound, buffer);
        return buffer.data();
    }

    detail::osqp_settings make_settings() const
    {
        detail::osqp_settings settings{OSQPSettings_new()};
        if(settings)
        {
            settings->eps_abs = m_eps_abs;
            settings->eps_rel = m_eps_rel;
            settings->max_iter = m_max_iter;
            settings->verbose = m_verbose;
            settings->warm_starting = m_warm_starting;
            settings->polishing = m_polishing;
        }
        return settings;
    }

    ctrlpp::expected<void, osqp_setup_error> configure_and_create_solver(const qp_problem<double>& problem)
    {
        auto p_csc = detail::make_csc(m_p_upper);
        auto a_csc = detail::make_csc(m_a);

        const auto settings = make_settings();
        if(!settings)
            return ctrlpp::unexpected(osqp_setup_error::settings_allocation_failed);

        OSQPSolver* solver = nullptr;
        auto exit_flag =
                osqp_setup(&solver, &p_csc, const_cast<OSQPFloat*>(problem.q.data()), &a_csc, m_osqp_lower.data(), m_osqp_upper.data(), m_n_constraints, m_n_vars, settings.get());
        m_solver.reset(solver);

        if(exit_flag != 0)
        {
            // osqp_setup publishes the solver pointer before the allocations that can fail, and its error paths free nothing,
            // so a failed setup hands back a partially built solver that is ours to release. Dropping it instead would leak it.
            m_solver.reset();
            return ctrlpp::unexpected(osqp_setup_error::setup_failed);
        }

        return {};
    }
};

}

#endif
