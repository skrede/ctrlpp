#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_PROBE_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_PROBE_H

#include "registry.h"
#include "hpipm_arm.h"
#include "json_writer.h"

#include <string>
#include <vector>
#include <optional>

namespace ctrlpp::bench::qp_solvers
{

struct probe_result
{
    std::string arm;
    bool passed;
    double max_dual_error;
    std::optional<double> solve_tolerance;
};

/// The HPIPM signs under which its probe answer matches the exact multipliers, or nothing when no
/// single sign pair does; error is the dual error under the returned (or best) pair.
struct hpipm_settlement
{
    std::optional<hpipm_dual_signs> signs;
    hpipm_dual_signs best;
    double error;
    std::int32_t matching_pairs;
};

auto settle_hpipm_signs() -> hpipm_settlement;

/// Each native-dual arm of the cohort solves a tiny problem whose exact primal-dual pair is known by
/// enumeration, at its high-tier settings (the mpc tier for the report-only arms; ProxQP at
/// proxqp_probe_solver_tolerance), and its mapped multipliers are compared with the exact ones.
auto run_probes(cohort group, const hpipm_settlement& hpipm) -> std::vector<probe_result>;

void write_probes(json_writer& out, const std::vector<probe_result>& results, const hpipm_settlement& hpipm);

}

#endif
