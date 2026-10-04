#include "hpipm_arm.h"
#include "hpipm_stages.h"
#include "hpipm_settings.h"

#include <array>
#include <vector>
#include <cstddef>
#include <string_view>

namespace ctrlpp::bench::qp_solvers
{

auto config_of(tier level) -> hpipm_config
{
    constexpr int mpc_iter_max = 50;
    constexpr int high_iter_max = 100;
    constexpr int warm_start_primal = 1;
    if(level == tier::mpc)
        return {SPEED, solver_tolerance(level), mpc_iter_max, warm_start_primal};
    return {BALANCE, solver_tolerance(level), high_iter_max, warm_start_primal};
}

auto reference_config() -> hpipm_config
{
    hpipm_config config = config_of(tier::high);
    config.tolerance = reference_solver_tolerance;
    return config;
}

void configure(const hpipm_config& config, d_ocp_qp_ipm_arg& arg)
{
    double tolerance = config.tolerance;
    int iter_max = config.iter_max;
    int warm_start = config.warm_start;
    d_ocp_qp_ipm_arg_set_default(config.mode, &arg);
    d_ocp_qp_ipm_arg_set_tol_stat(&tolerance, &arg);
    d_ocp_qp_ipm_arg_set_tol_eq(&tolerance, &arg);
    d_ocp_qp_ipm_arg_set_tol_ineq(&tolerance, &arg);
    d_ocp_qp_ipm_arg_set_tol_comp(&tolerance, &arg);
    d_ocp_qp_ipm_arg_set_iter_max(&iter_max, &arg);
    d_ocp_qp_ipm_arg_set_warm_start(&warm_start, &arg);
}

auto hpipm_status_token(int status) -> std::string_view
{
    constexpr std::array<std::string_view, 7> names{"SUCCESS", "MAX_ITER", "MIN_STEP", "NAN_SOL", "INCONS_EQ", "FAILED_REG",
                                                    "IN_PROGRESS"};
    return status >= 0 && status < static_cast<int>(names.size()) ? names[static_cast<std::size_t>(status)] : "UNKNOWN";
}

auto hpipm_claim(int status) -> solver_claim
{
    if(status == SUCCESS)
        return solver_claim::solved;
    if(status == MAX_ITER)
        return solver_claim::iteration_cap;
    return status == INCONS_EQ ? solver_claim::infeasible : solver_claim::error;
}

namespace
{

void write_defaults(json_writer& out, const d_ocp_qp_ipm_arg& arg)
{
    out.number("mu0", arg.mu0);
    out.number("alpha_min", arg.alpha_min);
    out.number("reg_prim", arg.reg_prim);
    out.integer("pred_corr", arg.pred_corr);
    out.integer("cond_pred_corr", arg.cond_pred_corr);
    out.integer("itref_pred_max", arg.itref_pred_max);
    out.integer("itref_corr_max", arg.itref_corr_max);
    out.integer("square_root_alg", arg.square_root_alg);
    out.integer("lq_fact", arg.lq_fact);
    out.integer("abs_form", arg.abs_form);
    out.integer("comp_dual_sol_eq", arg.comp_dual_sol_eq);
    out.integer("comp_res_exit", arg.comp_res_exit);
    out.integer("split_step", arg.split_step);
    out.integer("var_init_scheme", arg.var_init_scheme);
    out.integer("t_lam_min", arg.t_lam_min);
    out.integer("t0_init", arg.t0_init);
}

}

/// The defaults of a mode do not depend on the problem's size, so a one-stage problem's argument
/// structure reports them.
void describe_hpipm(json_writer& out, tier level)
{
    const hpipm_config config = config_of(level);
    d_ocp_qp_dim dim;
    std::vector<std::byte> dim_memory(d_ocp_qp_dim_memsize(1));
    d_ocp_qp_dim_create(1, &dim, dim_memory.data());
    d_ocp_qp_ipm_arg arg;
    std::vector<std::byte> arg_memory(d_ocp_qp_ipm_arg_memsize(&dim));
    d_ocp_qp_ipm_arg_create(&dim, &arg, arg_memory.data());
    configure(config, arg);
    out.text("mode", config.mode == SPEED ? "SPEED" : "BALANCE");
    out.number("tol_stat", arg.res_g_max);
    out.number("tol_eq", arg.res_b_max);
    out.number("tol_ineq", arg.res_d_max);
    out.number("tol_comp", arg.res_m_max);
    out.integer("iter_max", arg.iter_max);
    out.integer("warm_start", arg.warm_start);
    write_defaults(out, arg);
}

}
