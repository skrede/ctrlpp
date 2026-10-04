#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_DESCRIBE_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_DESCRIBE_H

#include <string>

namespace ctrlpp::bench::qp_solvers
{

/// What the harness itself runs with, for the run script to check against the pre-registration and
/// fill the manifest from: every arm's effective settings per tier, its update API, the expressibility
/// matrix, the harness's constants, and the long-horizon terminal weights' entries.
auto write_description(const std::string& path) -> bool;

}

#endif
