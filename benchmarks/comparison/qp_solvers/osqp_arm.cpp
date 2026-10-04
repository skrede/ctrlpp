#include "osqp_arm.h"

#include <osqp.h>

#include <array>
#include <vector>
#include <string_view>

namespace ctrlpp::bench::qp_solvers
{

namespace
{

/// A CSC copy in OSQP's index type, which is 64-bit in its default build.
struct csc_storage
{
    std::vector<OSQPInt> outer;
    std::vector<OSQPInt> inner;
    std::vector<OSQPFloat> values;
    OSQPCscMatrix matrix;
};

void store_csc(const Eigen::SparseMatrix<double>& source, csc_storage& out)
{
    const Eigen::SparseMatrix<double> compressed = source;
    out.outer.assign(compressed.outerIndexPtr(), compressed.outerIndexPtr() + compressed.cols() + 1);
    out.inner.assign(compressed.innerIndexPtr(), compressed.innerIndexPtr() + compressed.nonZeros());
    out.values.assign(compressed.valuePtr(), compressed.valuePtr() + compressed.nonZeros());
    OSQPCscMatrix_set_data(&out.matrix, compressed.rows(), compressed.cols(), compressed.nonZeros(), out.values.data(),
                           out.inner.data(), out.outer.data());
}

auto osqp_settings(tier level, bool forced_polish) -> OSQPSettings
{
    OSQPSettings settings;
    osqp_set_default_settings(&settings);
    settings.eps_abs = solver_tolerance(level);
    settings.eps_rel = solver_tolerance(level);
    settings.max_iter = iteration_cap;
    settings.warm_starting = 1;
    settings.polishing = level == tier::high || forced_polish ? 1 : 0;
    settings.adaptive_rho = OSQP_ADAPTIVE_RHO_UPDATE_ITERATIONS;
    settings.adaptive_rho_interval = osqp_adaptive_rho_interval;
    settings.verbose = 0;
    return settings;
}

auto status_token(OSQPInt status) -> std::string_view
{
    constexpr std::array<std::string_view, 11> names{
            "OSQP_SOLVED",          "OSQP_SOLVED_INACCURATE",  "OSQP_PRIMAL_INFEASIBLE", "OSQP_PRIMAL_INFEASIBLE_INACCURATE",
            "OSQP_DUAL_INFEASIBLE", "OSQP_DUAL_INFEASIBLE_INACCURATE", "OSQP_MAX_ITER_REACHED", "OSQP_TIME_LIMIT_REACHED",
            "OSQP_NON_CVX",         "OSQP_SIGINT",             "OSQP_UNSOLVED"};
    return status >= OSQP_SOLVED && status <= OSQP_UNSOLVED ? names[static_cast<std::size_t>(status - OSQP_SOLVED)] : "OSQP_UNKNOWN";
}

auto claim_of(OSQPInt status) -> solver_claim
{
    if(status == OSQP_SOLVED || status == OSQP_SOLVED_INACCURATE)
        return solver_claim::solved;
    if(status == OSQP_MAX_ITER_REACHED)
        return solver_claim::iteration_cap;
    if(status >= OSQP_PRIMAL_INFEASIBLE && status <= OSQP_DUAL_INFEASIBLE_INACCURATE)
        return solver_claim::infeasible;
    return solver_claim::error;
}

struct solver_cleanup
{
    void operator()(OSQPSolver* solver) const
    {
        osqp_cleanup(solver);
    }
};

class osqp_arm final : public qp_arm
{
public:
    explicit osqp_arm(const OSQPSettings& settings)
        : P_{}
        , A_{}
        , setup_error_(0)
        , settings_(settings)
        , solver_{}
    {
    }

    void setup(const qp_instance& instance, const qp_tick& tick) override
    {
        const Eigen::SparseMatrix<double> upper = instance.P_sparse.triangularView<Eigen::Upper>();
        store_csc(upper, P_);
        store_csc(instance.A_sparse, A_);
        OSQPSolver* solver = nullptr;
        setup_error_ = osqp_setup(&solver, &P_.matrix, tick.q.data(), &A_.matrix, tick.l.data(), tick.u.data(),
                                  instance.A_sparse.rows(), instance.A_sparse.cols(), &settings_);
        solver_.reset(solver);
    }

    void update(const qp_tick& tick) override
    {
        if(setup_error_ == 0)
            osqp_update_data_vec(solver_.get(), tick.q.data(), tick.l.data(), tick.u.data());
    }

    void solve() override
    {
        if(setup_error_ == 0)
            osqp_solve(solver_.get());
    }

    void answer(arm_answer& out) override;

private:
    csc_storage P_;
    csc_storage A_;
    OSQPInt setup_error_;
    OSQPSettings settings_;
    std::unique_ptr<OSQPSolver, solver_cleanup> solver_;
};

void osqp_arm::answer(arm_answer& out)
{
    out.fallback = false;
    if(setup_error_ != 0)
    {
        out = {{}, {}, -1, "OSQP_SETUP_ERROR", solver_claim::error, false};
        return;
    }
    OSQPInt m = 0;
    OSQPInt n = 0;
    osqp_get_dimensions(solver_.get(), &m, &n);
    out.x = Eigen::Map<const Eigen::VectorXd>(solver_->solution->x, n);
    out.y = Eigen::Map<const Eigen::VectorXd>(solver_->solution->y, m);
    out.iterations = static_cast<std::int32_t>(solver_->info->iter);
    out.status = status_token(solver_->info->status_val);
    out.claim = claim_of(solver_->info->status_val);
}

}

auto make_osqp_arm(tier level, bool forced_polish) -> std::unique_ptr<qp_arm>
{
    return std::make_unique<osqp_arm>(osqp_settings(level, forced_polish));
}

namespace
{

void write_defaults(json_writer& out, const OSQPSettings& s)
{
    out.text("adaptive_rho_rule", s.adaptive_rho == OSQP_ADAPTIVE_RHO_UPDATE_ITERATIONS ? "iterations" : "other");
    out.number("adaptive_rho_fraction", s.adaptive_rho_fraction);
    out.number("adaptive_rho_tolerance", s.adaptive_rho_tolerance);
    out.number("rho", s.rho);
    out.flag("rho_is_vec", s.rho_is_vec != 0);
    out.number("sigma", s.sigma);
    out.number("alpha", s.alpha);
    out.integer("scaling", s.scaling);
    out.integer("linsys_solver", s.linsys_solver);
    out.number("eps_prim_inf", s.eps_prim_inf);
    out.number("eps_dual_inf", s.eps_dual_inf);
    out.flag("scaled_termination", s.scaled_termination != 0);
    out.integer("check_termination", s.check_termination);
    out.flag("check_dualgap", s.check_dualgap != 0);
    out.number("time_limit", s.time_limit);
    out.number("delta", s.delta);
    out.integer("polish_refine_iter", s.polish_refine_iter);
}

}

void describe_osqp(json_writer& out, tier level, bool forced_polish)
{
    const OSQPSettings s = osqp_settings(level, forced_polish);
    out.number("eps_abs", s.eps_abs);
    out.number("eps_rel", s.eps_rel);
    out.integer("max_iter", s.max_iter);
    out.flag("warm_starting", s.warm_starting != 0);
    out.flag("polishing", s.polishing != 0);
    out.flag("adaptive_rho", s.adaptive_rho != OSQP_ADAPTIVE_RHO_UPDATE_DISABLED);
    out.integer("adaptive_rho_interval", s.adaptive_rho_interval);
    out.flag("verbose", s.verbose != 0);
    write_defaults(out, s);
}

}
