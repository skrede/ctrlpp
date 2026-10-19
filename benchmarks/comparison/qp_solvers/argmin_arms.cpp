#include "argmin_arms.h"
#include "argmin_outcome.h"

#include <argmin/qp/dense_admm_qp.h>
#include <argmin/qp/sparse_admm_qp.h>
#include <argmin/qp/dense_active_set_qp.h>

#include <optional>

namespace ctrlpp::bench::qp_solvers
{

namespace
{

/// The solve call takes the vectors itself, so update() stages them into buffers the solver reads.
template <typename Solver, typename Options>
class argmin_arm_base : public qp_arm
{
public:
    explicit argmin_arm_base(const Options& options)
        : options_(options)
        , q_{}
        , l_{}
        , u_{}
        , solver_{}
        , result_{}
        , error_{}
    {
    }

    void update(const qp_tick& tick) override
    {
        stage(tick);
    }

    void solve() override
    {
        error_ = solver_->resolve_into(q_, l_, u_, result_, options_);
    }

    void answer(arm_answer& out) override
    {
        fill_answer(result_, error_, out);
        out.fallback = false;
    }

protected:
    Options options_;
    Eigen::VectorXd q_;
    Eigen::VectorXd l_;
    Eigen::VectorXd u_;
    std::optional<Solver> solver_;
    typename Solver::result_type result_;
    std::optional<argmin::qp_error> error_;

    void stage(const qp_tick& tick)
    {
        q_ = tick.q;
        l_ = tick.l;
        u_ = tick.u;
    }
};

using dense_admm = argmin::dense_admm_qp_solver<double, Eigen::Dynamic>;
using dense_active_set = argmin::dense_active_set_qp_solver<double, Eigen::Dynamic>;
using sparse_admm = argmin::sparse_admm_qp_solver<double>;

class dense_admm_arm final : public argmin_arm_base<dense_admm, argmin::dense_qp_options<double>>
{
public:
    using argmin_arm_base::argmin_arm_base;

    void setup(const qp_instance& instance, const qp_tick& tick) override
    {
        stage(tick);
        solver_.emplace(static_cast<int>(instance.P.rows()), static_cast<int>(instance.A.rows()));
        error_ = solver_->solve_into(instance.P, q_, instance.A, l_, u_, result_, options_);
    }
};

class sparse_admm_arm final : public argmin_arm_base<sparse_admm, argmin::sparse_qp_options<double>>
{
public:
    using argmin_arm_base::argmin_arm_base;

    void setup(const qp_instance& instance, const qp_tick& tick) override
    {
        stage(tick);
        solver_.emplace();
        error_ = solver_->solve_into(instance.P_sparse, q_, instance.A_sparse, l_, u_, result_, options_);
    }
};

/// The active-set method needs a feasible start; when the retained solution is infeasible for the new
/// bounds the solver is reset to U = 0 (feasible for every posed tick) and the resolve repeated.
class active_set_arm final : public argmin_arm_base<dense_active_set, argmin::dense_active_set_qp_options<double>>
{
public:
    explicit active_set_arm(const argmin::dense_active_set_qp_options<double>& options)
        : argmin_arm_base(options)
        , fallback_(false)
    {
    }

    void setup(const qp_instance& instance, const qp_tick& tick) override
    {
        stage(tick);
        solver_.emplace(static_cast<int>(instance.P.rows()), static_cast<int>(instance.A.rows()));
        error_ = solver_->solve_into(instance.P, q_, instance.A, l_, u_, result_, options_);
    }

    void solve() override
    {
        fallback_ = false;
        error_ = solver_->resolve_into(q_, l_, u_, result_, options_);
        if(error_ != argmin::qp_error::infeasible_start)
            return;
        fallback_ = true;
        solver_->reset();
        error_ = solver_->resolve_into(q_, l_, u_, result_, options_);
    }

    void answer(arm_answer& out) override
    {
        fill_answer(result_, error_, out);
        out.fallback = fallback_;
    }

private:
    bool fallback_;
};

auto active_set_options(double tolerance) -> argmin::dense_active_set_qp_options<double>
{
    argmin::dense_active_set_qp_options<double> options;
    options.tolerance = tolerance;
    options.max_iterations = iteration_cap;
    return options;
}

}

auto make_argmin_dense_admm_arm(tier level, bool forced_polish) -> std::unique_ptr<qp_arm>
{
    return std::make_unique<dense_admm_arm>(admm_options<argmin::dense_qp_options<double>>(level, forced_polish));
}

auto make_argmin_active_set_arm(tier level) -> std::unique_ptr<qp_arm>
{
    return std::make_unique<active_set_arm>(active_set_options(solver_tolerance(level)));
}

auto make_argmin_sparse_admm_arm(tier level, bool forced_polish) -> std::unique_ptr<qp_arm>
{
    return std::make_unique<sparse_admm_arm>(admm_options<argmin::sparse_qp_options<double>>(level, forced_polish));
}

auto make_argmin_active_set_reference() -> std::unique_ptr<qp_arm>
{
    return std::make_unique<active_set_arm>(active_set_options(reference_solver_tolerance));
}

}
