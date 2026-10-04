#include "argmin_arms.h"
#include "argmin_outcome.h"

#include <argmin/options/dense_qp_options.h>
#include <argmin/options/sparse_qp_options.h>
#include <argmin/options/dense_active_set_qp_options.h>

namespace ctrlpp::bench::qp_solvers
{

namespace
{

template <typename Options>
void write_admm_options(json_writer& out, const Options& options)
{
    out.number("eps_abs", options.eps_abs);
    out.number("eps_rel", options.eps_rel);
    out.integer("max_iterations", options.max_iterations);
    out.flag("warm_start", options.warm_start);
    out.flag("polish", options.polish);
    out.number("rho", options.rho);
    out.number("sigma", options.sigma);
    out.number("alpha", options.alpha);
    out.number("eps_prim_inf", options.eps_prim_inf);
    out.number("eps_dual_inf", options.eps_dual_inf);
    out.integer("check_termination", options.check_termination);
    out.flag("adaptive_rho", options.adaptive_rho);
    out.integer("adaptive_rho_interval", options.adaptive_rho_interval);
    out.number("adaptive_rho_tolerance", options.adaptive_rho_tolerance);
    out.integer("scaling", options.scaling);
    out.integer("polish_refine_iter", options.polish_refine_iter);
    out.number("delta", options.delta);
    out.flag("validate_convexity", options.validate_convexity);
}

}

void describe_argmin_dense_admm(json_writer& out, tier level, bool forced_polish)
{
    write_admm_options(out, admm_options<argmin::dense_qp_options>(level, forced_polish));
}

void describe_argmin_sparse_admm(json_writer& out, tier level, bool forced_polish)
{
    write_admm_options(out, admm_options<argmin::sparse_qp_options>(level, forced_polish));
}

void describe_argmin_active_set(json_writer& out, tier level)
{
    argmin::dense_active_set_qp_options options;
    options.tolerance = solver_tolerance(level);
    options.max_iterations = iteration_cap;
    out.number("tolerance", options.tolerance);
    out.integer("max_iterations", options.max_iterations);
    out.number("equality_tolerance", options.equality_tolerance);
}

}
