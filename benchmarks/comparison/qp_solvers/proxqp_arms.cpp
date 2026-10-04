#include "row_split.h"
#include "competitor_arms.h"

#include <proxsuite/proxqp/dense/dense.hpp>
#include <proxsuite/proxqp/sparse/sparse.hpp>

#include <array>
#include <optional>
#include <string_view>

namespace ctrlpp::bench::qp_solvers
{

namespace
{

namespace proxqp = proxsuite::proxqp;

template <typename Settings>
void apply_tier(Settings& settings, tier level)
{
    settings.eps_abs = solver_tolerance(level);
    settings.eps_rel = solver_tolerance(level);
    settings.max_iter = iteration_cap;
}

auto status_token(proxqp::QPSolverOutput status) -> std::string_view
{
    constexpr std::array<std::string_view, 6> names{"PROXQP_SOLVED",
                                                    "PROXQP_MAX_ITER_REACHED",
                                                    "PROXQP_PRIMAL_INFEASIBLE",
                                                    "PROXQP_SOLVED_CLOSEST_PRIMAL_FEASIBLE",
                                                    "PROXQP_DUAL_INFEASIBLE",
                                                    "PROXQP_NOT_RUN"};
    return names[static_cast<std::size_t>(status)];
}

auto claim_of(proxqp::QPSolverOutput status) -> solver_claim
{
    switch(status)
    {
        case proxqp::QPSolverOutput::PROXQP_SOLVED:
        case proxqp::QPSolverOutput::PROXQP_SOLVED_CLOSEST_PRIMAL_FEASIBLE:
            return solver_claim::solved;
        case proxqp::QPSolverOutput::PROXQP_MAX_ITER_REACHED:
            return solver_claim::iteration_cap;
        case proxqp::QPSolverOutput::PROXQP_PRIMAL_INFEASIBLE:
        case proxqp::QPSolverOutput::PROXQP_DUAL_INFEASIBLE:
            return solver_claim::infeasible;
        default:
            return solver_claim::error;
    }
}

/// Dense and sparse ProxQP share one wrapper interface; only the matrix types handed to init() differ.
/// A solve that warm-starts from a previous result before any result exists returns NaN, so setup()
/// solves once with ProxQP's default initial guess and every later solve warm-starts.
template <typename Solver>
class proxqp_arm final : public qp_arm
{
public:
    explicit proxqp_arm(tier level)
        : level_(level)
        , rows_{}
        , b_{}
        , lower_{}
        , upper_{}
        , solver_{}
    {
    }

    void setup(const qp_instance& instance, const qp_tick& tick) override
    {
        rows_ = split_rows(tick.l, tick.u);
        const Eigen::Index n = instance.P.rows();
        solver_.emplace(n, rows_.equalities.size(), rows_.inequalities.size());
        apply_tier(solver_->settings, level_);
        configure();
        rows_.gather(tick.l, tick.u, b_, lower_, upper_);
        initialize(instance, tick);
        solver_->solve();
        solver_->settings.initial_guess = proxqp::InitialGuessStatus::WARM_START_WITH_PREVIOUS_RESULT;
    }

    void update(const qp_tick& tick) override
    {
        rows_.gather(tick.l, tick.u, b_, lower_, upper_);
        solver_->update(std::nullopt, tick.q, std::nullopt, b_, std::nullopt, lower_, upper_);
    }

    void solve() override
    {
        solver_->solve();
    }

    void answer(arm_answer& out) override
    {
        out.x = solver_->results.x;
        rows_.scatter(solver_->results.y, solver_->results.z, out.y);
        out.iterations = static_cast<std::int32_t>(solver_->results.info.iter);
        out.status = status_token(solver_->results.info.status);
        out.claim = claim_of(solver_->results.info.status);
        out.fallback = false;
    }

private:
    tier level_;
    row_split rows_;
    Eigen::VectorXd b_;
    Eigen::VectorXd lower_;
    Eigen::VectorXd upper_;
    std::optional<Solver> solver_;

    void configure();

    void initialize(const qp_instance& instance, const qp_tick& tick);
};

using dense_qp = proxqp::dense::QP<double>;
using sparse_qp = proxqp::sparse::QP<double, int>;

template <>
void proxqp_arm<dense_qp>::configure()
{
}

template <>
void proxqp_arm<sparse_qp>::configure()
{
    solver_->settings.sparse_backend = proxqp::SparseBackend::SparseCholesky;
}

template <>
void proxqp_arm<dense_qp>::initialize(const qp_instance& instance, const qp_tick& tick)
{
    const Eigen::MatrixXd A = rows_.equality_rows(instance.A);
    const Eigen::MatrixXd C = rows_.inequality_rows(instance.A);
    solver_->init(instance.P, tick.q, A, b_, C, lower_, upper_);
}

template <>
void proxqp_arm<sparse_qp>::initialize(const qp_instance& instance, const qp_tick& tick)
{
    const Eigen::SparseMatrix<double> A = rows_.equality_rows(instance.A).sparseView();
    const Eigen::SparseMatrix<double> C = rows_.inequality_rows(instance.A).sparseView();
    solver_->init(instance.P_sparse, tick.q, A, b_, C, lower_, upper_);
}

}

auto make_proxqp_dense_arm(tier level) -> std::unique_ptr<qp_arm>
{
    return std::make_unique<proxqp_arm<dense_qp>>(level);
}

auto make_proxqp_sparse_arm(tier level) -> std::unique_ptr<qp_arm>
{
    return std::make_unique<proxqp_arm<sparse_qp>>(level);
}

void describe_proxqp(json_writer& out, tier level, bool sparse)
{
    proxqp::Settings<double> s;
    apply_tier(s, level);
    out.number("eps_abs", s.eps_abs);
    out.number("eps_rel", s.eps_rel);
    out.integer("max_iter", s.max_iter);
    out.text("initial_guess", "WARM_START_WITH_PREVIOUS_RESULT");
    out.text("setup_solve_initial_guess", "EQUALITY_CONSTRAINED_INITIAL_GUESS");
    if(sparse)
        out.text("backend", "SparseCholesky");
    out.number("default_rho", s.default_rho);
    out.number("default_mu_eq", s.default_mu_eq);
    out.number("default_mu_in", s.default_mu_in);
    out.number("alpha_bcl", s.alpha_bcl);
    out.number("beta_bcl", s.beta_bcl);
    out.number("eps_primal_inf", s.eps_primal_inf);
    out.number("eps_dual_inf", s.eps_dual_inf);
    out.integer("nb_iterative_refinement", s.nb_iterative_refinement);
    out.integer("preconditioner_max_iter", s.preconditioner_max_iter);
    out.flag("compute_preconditioner", s.compute_preconditioner);
    out.flag("compute_timings", s.compute_timings);
    out.flag("check_duality_gap", s.check_duality_gap);
}

}
