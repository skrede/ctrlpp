#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_TINYMPC_ARM_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_TINYMPC_ARM_H

#include "arm.h"
#include "registry.h"
#include "json_writer.h"

#include <memory>
#include <string>

namespace ctrlpp::bench::qp_solvers
{

/// TinyMPC exposes no usable multipliers, so its answer carries none; the harness recovers them.
auto make_tinympc_arm(tier level) -> std::unique_ptr<qp_arm>;

void describe_tinympc_arm(json_writer& out, tier level);

auto tinympc_unit_eigen_version() -> std::string;

}

#endif
