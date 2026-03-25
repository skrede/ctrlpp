#ifndef HPP_GUARD_CTRLPP_MPC_OSQP_SOLVER_H
#define HPP_GUARD_CTRLPP_MPC_OSQP_SOLVER_H

#include "ctrlpp/mpc/qp_types.h"

#include <Eigen/Sparse>

#include <osqp.h>

#include <stdexcept>
#include <type_traits>
#include <utility>

namespace ctrlpp
{

class osqp_solver
{
    static_assert(std::is_same_v<double, OSQPFloat>, "OSQP must be compiled with double precision");

    using sparse_t = Eigen::SparseMatrix<OSQPFloat, Eigen::ColMajor, OSQPInt>;

public:
    using scalar_type = double;

    explicit osqp_solver(double eps_abs = 1e-3, double eps_rel = 1e-3, int max_iter = 4000, bool verbose = false, bool warm_starting = true, bool polishing = true)
        : eps_abs_{eps_abs}, eps_rel_{eps_rel}, max_iter_{max_iter}, verbose_{verbose}, warm_starting_{warm_starting}, polishing_{polishing}
    {
    }

    ~osqp_solver() { cleanup(); }

    osqp_solver(const osqp_solver&) = delete;
    auto operator=(const osqp_solver&) -> osqp_solver& = delete;

    osqp_solver(osqp_solver&& other) noexcept
        : solver_{std::exchange(other.solver_, nullptr)}
        , p_upper_{std::move(other.p_upper_)}
        , a_{std::move(other.a_)}
        , n_{other.n_}
        , m_{other.m_}
        , eps_abs_{other.eps_abs_}
        , eps_rel_{other.eps_rel_}
        , max_iter_{other.max_iter_}
        , verbose_{other.verbose_}
        , warm_starting_{other.warm_starting_}
        , polishing_{other.polishing_}
    {
    }

    auto operator=(osqp_solver&& other) noexcept -> osqp_solver&
    {
        if(this != &other)
        {
            cleanup();
            solver_ = std::exchange(other.solver_, nullptr);
            p_upper_ = std::move(other.p_upper_);
            a_ = std::move(other.a_);
            n_ = other.n_;
            m_ = other.m_;
            eps_abs_ = other.eps_abs_;
            eps_rel_ = other.eps_rel_;
            max_iter_ = other.max_iter_;
            verbose_ = other.verbose_;
            warm_starting_ = other.warm_starting_;
            polishing_ = other.polishing_;
        }
        return *this;
    }

    void setup(const qp_problem<double>& problem)
    {
        cleanup();
        prepare_sparse_matrices(problem);
        configure_and_create_solver(problem);
    }

    auto solve(const qp_update<double>& update) -> qp_result<double>
    {
        osqp_update_data_vec(solver_, update.q.size() > 0 ? update.q.data() : nullptr, update.l.size() > 0 ? update.l.data() : nullptr, update.u.size() > 0 ? update.u.data() : nullptr);

        if(update.warm_x.size() > 0 || update.warm_y.size() > 0)
        {
            osqp_warm_start(solver_, update.warm_x.size() > 0 ? update.warm_x.data() : nullptr, update.warm_y.size() > 0 ? update.warm_y.data() : nullptr);
        }

        osqp_solve(solver_);

        qp_result<double> result;
        result.x = Eigen::Map<const Eigen::VectorXd>(solver_->solution->x, n_);
        result.y = Eigen::Map<const Eigen::VectorXd>(solver_->solution->y, m_);
        result.status = map_status(solver_->info->status_val);
        result.objective = solver_->info->obj_val;
        result.solve_time = solver_->info->solve_time;
        result.iterations = static_cast<int>(solver_->info->iter);
        result.primal_residual = solver_->info->prim_res;
        result.dual_residual = solver_->info->dual_res;

        return result;
    }

private:
    void prepare_sparse_matrices(const qp_problem<double>& problem)
    {
        p_upper_ = problem.P.template triangularView<Eigen::Upper>();
        p_upper_.makeCompressed();

        a_ = problem.A.template cast<OSQPFloat>();
        a_.makeCompressed();

        n_ = static_cast<OSQPInt>(problem.P.cols());
        m_ = static_cast<OSQPInt>(problem.A.rows());
    }

    void configure_and_create_solver(const qp_problem<double>& problem)
    {
        auto p_csc = make_csc(p_upper_);
        auto a_csc = make_csc(a_);

        auto* settings = OSQPSettings_new();
        if(!settings)
            throw std::runtime_error("OSQP settings allocation failed");

        settings->eps_abs = eps_abs_;
        settings->eps_rel = eps_rel_;
        settings->max_iter = max_iter_;
        settings->verbose = verbose_;
        settings->warm_starting = warm_starting_;
        settings->polishing = polishing_;

        auto exit_flag = osqp_setup(&solver_, &p_csc, const_cast<OSQPFloat*>(problem.q.data()), &a_csc, const_cast<OSQPFloat*>(problem.l.data()), const_cast<OSQPFloat*>(problem.u.data()), m_, n_, settings);

        OSQPSettings_free(settings);

        if(exit_flag != 0)
        {
            solver_ = nullptr;
            throw std::runtime_error("OSQP setup failed");
        }
    }

    void cleanup()
    {
        if(solver_)
        {
            osqp_cleanup(solver_);
            solver_ = nullptr;
        }
    }

    static auto make_csc(const sparse_t& mat) -> OSQPCscMatrix
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

    static constexpr auto map_status(OSQPInt status_val) -> solve_status
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

    OSQPSolver* solver_{nullptr};
    sparse_t p_upper_;
    sparse_t a_;
    OSQPInt n_{0};
    OSQPInt m_{0};

    double eps_abs_;
    double eps_rel_;
    int max_iter_;
    bool verbose_;
    bool warm_starting_;
    bool polishing_;
};

} // namespace ctrlpp

#endif
