#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_HPIPM_SETTINGS_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_HPIPM_SETTINGS_H

#include "arm.h"
#include "registry.h"

extern "C"
{
#include <hpipm_common.h>
#include <hpipm_d_ocp_qp_sol.h>
#include <hpipm_d_ocp_qp_ipm.h>
}

#include <string_view>

namespace ctrlpp::bench::qp_solvers
{

/// One tolerance for the stationarity, equality, inequality and complementarity residuals.
struct hpipm_config
{
    hpipm_mode mode;
    double tolerance;
    int iter_max;
    int warm_start;
};

auto config_of(tier level) -> hpipm_config;

auto reference_config() -> hpipm_config;

void configure(const hpipm_config& config, d_ocp_qp_ipm_arg& arg);

auto hpipm_status_token(int status) -> std::string_view;

auto hpipm_claim(int status) -> solver_claim;

}

#endif
