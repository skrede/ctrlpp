#include "json_writer.h"
#include "tinympc_session.h"

#include <tinympc/tiny_api.hpp>

#include <format>
#include <limits>
#include <iostream>

namespace ctrlpp::bench::qp_solvers
{

namespace
{

auto as_matrix(const std::vector<double>& values, std::int32_t rows, std::int32_t cols) -> tinyMatrix
{
    return Eigen::Map<const tinyMatrix>(values.data(), rows, cols);
}

/// Unbounded except, with a state box, the actuated states of knots 1..N.
auto state_bounds(const tinympc_problem& p, double sign) -> tinyMatrix
{
    const double infinity = std::numeric_limits<double>::infinity();
    tinyMatrix bound = tinyMatrix::Constant(p.nx, p.horizon + 1, sign * infinity);
    if(p.state_limit)
        bound.block(0, 1, p.nu, p.horizon).setConstant(sign * *p.state_limit);
    return bound;
}

void apply_options(const tinympc_options& options, bool state_box, TinySettings& settings)
{
    settings.abs_pri_tol = options.tolerance;
    settings.abs_dua_tol = options.tolerance;
    settings.max_iter = options.max_iter;
    settings.adaptive_rho = options.adaptive_rho ? 1 : 0;
    settings.en_state_bound = state_box ? 1 : 0;
    settings.en_input_bound = 1;
}

/// tiny_setup allocates its parts with new and TinyMPC offers no call to release them.
struct tiny_release
{
    void operator()(TinySolver* solver) const
    {
        delete solver->solution;
        delete solver->settings;
        delete solver->cache;
        delete solver->work;
        delete solver;
    }
};

}

struct tinympc_session::solver_state
{
    tinyVector x0;
    std::unique_ptr<TinySolver, tiny_release> solver;
};

/// TinyMPC writes a line to std::cout on every converged solve. The stream is put in a failed state, so
/// that write formats and flushes nothing inside the timed solve; the harness itself never uses it.
tinympc_session::tinympc_session(const tinympc_problem& p, const tinympc_options& options)
    : state_(std::make_unique<solver_state>())
{
    std::cout.setstate(std::ios_base::badbit);
    const tinyMatrix A = as_matrix(p.transition, p.nx, p.nx);
    const tinyMatrix B = as_matrix(p.actuation, p.nx, p.nu);
    const tinyMatrix Q = Eigen::Map<const tinyVector>(p.state_weight.data(), p.nx).asDiagonal();
    const tinyMatrix R = Eigen::Map<const tinyVector>(p.input_weight.data(), p.nu).asDiagonal();
    const tinyMatrix f = tinyMatrix::Zero(p.nx, 1);
    TinySolver* solver = nullptr;
    tiny_setup(&solver, A, B, f, Q, R, options.rho, p.nx, p.nu, p.horizon + 1, 0);
    state_->solver.reset(solver);
    apply_options(options, p.state_limit.has_value(), *solver->settings);
    const tinyMatrix u_bound = tinyMatrix::Constant(p.nu, p.horizon, p.input_limit);
    tiny_set_bound_constraints(solver, state_bounds(p, -1.0), state_bounds(p, 1.0), -u_bound, u_bound);
    if(p.state_limit)
        tiny_set_x_ref(solver, as_matrix(p.target, p.nx, 1).replicate(1, p.horizon + 1));
    state_->x0 = tinyVector::Zero(p.nx);
}

tinympc_session::~tinympc_session() = default;

void tinympc_session::set_initial_state(const double* x0)
{
    state_->x0 = Eigen::Map<const tinyVector>(x0, state_->x0.size());
    tiny_set_x0(state_->solver.get(), state_->x0);
}

auto tinympc_session::solve() -> std::int32_t
{
    return tiny_solve(state_->solver.get());
}

auto tinympc_session::iterations() const -> std::int32_t
{
    return state_->solver->solution->iter;
}

void tinympc_session::read(double* states, double* inputs) const
{
    const TinySolution& solution = *state_->solver->solution;
    Eigen::Map<tinyMatrix>(states, solution.x.rows(), solution.x.cols()) = solution.x;
    Eigen::Map<tinyMatrix>(inputs, solution.u.rows(), solution.u.cols()) = solution.u;
}

auto tinympc_eigen_version() -> std::string
{
    return std::format("{}.{}.{}", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION);
}

/// The state-bound switch follows the variant rather than the tier, so it is not listed.
void describe_tinympc(json_writer& out, const tinympc_options& options)
{
    TinySettings s;
    tiny_set_default_settings(&s);
    apply_options(options, false, s);
    out.number("abs_pri_tol", s.abs_pri_tol);
    out.number("abs_dua_tol", s.abs_dua_tol);
    out.integer("max_iter", s.max_iter);
    out.number("rho", options.rho);
    out.flag("adaptive_rho", s.adaptive_rho != 0);
    out.integer("check_termination", s.check_termination);
    out.flag("en_input_bound", s.en_input_bound != 0);
    out.number("adaptive_rho_min", s.adaptive_rho_min);
    out.number("adaptive_rho_max", s.adaptive_rho_max);
    out.flag("adaptive_rho_enable_clipping", s.adaptive_rho_enable_clipping != 0);
}

}
