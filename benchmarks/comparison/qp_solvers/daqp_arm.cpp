#include "competitor_arms.h"

extern "C"
{
#include <api.h>
#include <utils.h>
}

#include <array>
#include <vector>
#include <string_view>

namespace ctrlpp::bench::qp_solvers
{

namespace
{

using row_major = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

auto daqp_settings(tier level) -> DAQPSettings
{
    DAQPSettings settings;
    daqp_default_settings(&settings);
    settings.primal_tol = solver_tolerance(level);
    settings.iter_limit = iteration_cap;
    return settings;
}

auto status_token(int exitflag) -> std::string_view
{
    constexpr std::array<std::string_view, 13> names{"DAQP_EXIT_UNSUPPORTED", "DAQP_EXIT_TIMELIMIT", "DAQP_EXIT_OVERDETERMINED_INITIAL",
                                                     "DAQP_EXIT_NONCONVEX",   "DAQP_EXIT_ITERLIMIT", "DAQP_EXIT_UNBOUNDED",
                                                     "DAQP_EXIT_CYCLE",       "DAQP_EXIT_INFEASIBLE", "DAQP_EXIT_UNKNOWN",
                                                     "DAQP_EXIT_OPTIMAL",     "DAQP_EXIT_SOFT_OPTIMAL", "DAQP_EXIT_UNKNOWN",
                                                     "DAQP_EXIT_OPTIMAL_INEXACT"};
    const int index = exitflag - DAQP_EXIT_UNSUPPORTED;
    return index >= 0 && index < static_cast<int>(names.size()) ? names[static_cast<std::size_t>(index)] : "DAQP_EXIT_UNKNOWN";
}

auto claim_of(int exitflag) -> solver_claim
{
    if(exitflag == DAQP_EXIT_OPTIMAL || exitflag == DAQP_EXIT_SOFT_OPTIMAL || exitflag == DAQP_EXIT_OPTIMAL_INEXACT)
        return solver_claim::solved;
    if(exitflag == DAQP_EXIT_ITERLIMIT)
        return solver_claim::iteration_cap;
    if(exitflag == DAQP_EXIT_INFEASIBLE || exitflag == DAQP_EXIT_UNBOUNDED)
        return solver_claim::infeasible;
    return solver_claim::error;
}

/// The workspace keeps its working set between solves; a new tick rewrites the problem's vectors in
/// place and refreshes the parts of the LDP that depend on them.
class daqp_arm final : public qp_arm
{
public:
    explicit daqp_arm(const DAQPSettings& settings)
        : H_{}
        , A_{}
        , setup_flag_(0)
        , f_{}
        , x_{}
        , result_{}
        , work_{}
        , problem_{}
        , settings_(settings)
        , upper_{}
        , lower_{}
        , lambda_{}
        , sense_{}
    {
    }

    ~daqp_arm() override
    {
        work_.settings = nullptr;
        free_daqp_workspace(&work_);
        free_daqp_ldp(&work_);
    }

    daqp_arm(const daqp_arm&) = delete;
    daqp_arm& operator=(const daqp_arm&) = delete;

    void setup(const qp_instance& instance, const qp_tick& tick) override;

    void update(const qp_tick& tick) override
    {
        f_ = tick.q;
        upper_ = tick.u;
        lower_ = tick.l;
        if(setup_flag_ >= 0)
            daqp_update_ldp(DAQP_UPDATE_v + DAQP_UPDATE_d, &work_, &problem_);
    }

    void solve() override
    {
        if(setup_flag_ >= 0)
            daqp_solve(&result_, &work_);
    }

    void answer(arm_answer& out) override
    {
        out.x = x_;
        out.y = lambda_;
        out.iterations = setup_flag_ >= 0 ? result_.iter : -1;
        out.status = setup_flag_ >= 0 ? status_token(result_.exitflag) : "DAQP_SETUP_ERROR";
        out.claim = setup_flag_ >= 0 ? claim_of(result_.exitflag) : solver_claim::error;
        out.fallback = false;
    }

private:
    row_major H_;
    row_major A_;
    int setup_flag_;
    Eigen::VectorXd f_;
    Eigen::VectorXd x_;
    DAQPResult result_;
    DAQPWorkspace work_;
    DAQPProblem problem_;
    DAQPSettings settings_;
    Eigen::VectorXd upper_;
    Eigen::VectorXd lower_;
    Eigen::VectorXd lambda_;
    std::vector<int> sense_;

    void pose(const qp_instance& instance);
};

void daqp_arm::pose(const qp_instance& instance)
{
    H_ = instance.P;
    A_ = instance.A.bottomRows(instance.A.rows() - instance.box_rows);
    problem_.n = static_cast<int>(instance.P.rows());
    problem_.m = static_cast<int>(instance.A.rows());
    problem_.ms = instance.box_rows;
    problem_.H = H_.data();
    problem_.f = f_.data();
    problem_.A = A_.data();
    problem_.bupper = upper_.data();
    problem_.blower = lower_.data();
    problem_.sense = sense_.data();
    result_.x = x_.data();
    result_.lam = lambda_.data();
}

void daqp_arm::setup(const qp_instance& instance, const qp_tick& tick)
{
    f_ = tick.q;
    upper_ = tick.u;
    lower_ = tick.l;
    x_.resize(instance.P.rows());
    lambda_.resize(instance.A.rows());
    sense_.assign(static_cast<std::size_t>(instance.A.rows()), 0);
    for(Eigen::Index i = 0; i < tick.l.size(); ++i)
        if(tick.l(i) == tick.u(i))
            sense_[static_cast<std::size_t>(i)] = DAQP_ACTIVE + DAQP_IMMUTABLE;
    pose(instance);
    work_.settings = &settings_;
    setup_flag_ = setup_daqp(&problem_, &work_, nullptr);
}

}

auto make_daqp_arm(tier level) -> std::unique_ptr<qp_arm>
{
    return std::make_unique<daqp_arm>(daqp_settings(level));
}

void describe_daqp(json_writer& out, tier level)
{
    const DAQPSettings s = daqp_settings(level);
    out.number("primal_tol", s.primal_tol);
    out.integer("iter_limit", s.iter_limit);
    out.number("dual_tol", s.dual_tol);
    out.number("zero_tol", s.zero_tol);
    out.number("pivot_tol", s.pivot_tol);
    out.number("progress_tol", s.progress_tol);
    out.integer("cycle_tol", s.cycle_tol);
    out.number("fval_bound", s.fval_bound);
    out.number("eps_prox", s.eps_prox);
    out.number("eta_prox", s.eta_prox);
    out.number("rho_soft", s.rho_soft);
    out.number("rel_subopt", s.rel_subopt);
    out.number("abs_subopt", s.abs_subopt);
    out.number("sing_tol", s.sing_tol);
    out.number("refactor_tol", s.refactor_tol);
    out.number("time_limit", s.time_limit);
    out.integer("eq_reduction", s.eq_reduction);
}

}
