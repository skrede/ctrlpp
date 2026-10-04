#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_OSQP_ARM_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_OSQP_ARM_H

#include "arm.h"
#include "registry.h"
#include "json_writer.h"

#include <memory>

namespace ctrlpp::bench::qp_solvers
{

/// Drives OSQP's C API directly rather than ctrlpp's wrapper, so the adaptive-rho interval can be pinned
/// to an iteration count instead of the wrapper's timing-based default.
auto make_osqp_arm(tier level, bool forced_polish) -> std::unique_ptr<qp_arm>;

void describe_osqp(json_writer& out, tier level, bool forced_polish);

}

#endif
