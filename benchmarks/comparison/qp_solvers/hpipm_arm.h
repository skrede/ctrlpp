#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_HPIPM_ARM_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_HPIPM_ARM_H

#include "arm.h"
#include "registry.h"
#include "json_writer.h"

#include <memory>

namespace ctrlpp::bench::qp_solvers
{

/// The signs that carry HPIPM's dynamics multipliers and its stage-0 state-bound multipliers into
/// OSQP's convention for the dynamics rows and the initial-state rows. They are settled by the
/// convention probe, not read from HPIPM's documentation.
struct hpipm_dual_signs
{
    double initial_state;
    double dynamics;
};

/// Poses the instance's optimal-control structure as an OCP QP with the initial state as an equality
/// bound on stage 0; each tick rewrites that bound.
auto make_hpipm_arm(tier level, hpipm_dual_signs signs) -> std::unique_ptr<qp_arm>;

auto make_hpipm_reference(hpipm_dual_signs signs) -> std::unique_ptr<qp_arm>;

void describe_hpipm(json_writer& out, tier level);

}

#endif
