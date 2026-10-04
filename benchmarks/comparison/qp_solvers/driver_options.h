#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_DRIVER_OPTIONS_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_DRIVER_OPTIONS_H

#include "registry.h"

#include <string>
#include <vector>
#include <cstdint>
#include <optional>
#include <string_view>

namespace ctrlpp::bench::qp_solvers
{

/// Each list filter is a comma-separated set of names; an empty filter admits everything.
/// on_probe_failure is "abort" unless a diagnostic smoke asks to "continue" past a failed convention
/// probe, which leaves the failure in the probe record for the verdict command to refuse.
struct driver_options
{
    cohort group;
    std::string out_dir;
    std::string run_id;
    std::int32_t rep;
    std::int32_t steps;
    std::vector<std::string> cells;
    std::vector<std::string> variants;
    std::vector<std::string> regimes;
    std::vector<std::string> tiers;
    std::vector<std::string> arms;
    std::vector<std::string> refused_families;
    std::string describe_path;
    std::string gate_path;
    std::string on_probe_failure;
};

/// Nothing when the arguments are malformed; the reason has then been printed.
auto parse_options(cohort group, int argc, char** argv) -> std::optional<driver_options>;

auto admits(const std::vector<std::string>& filter, std::string_view name) -> bool;

}

#endif
