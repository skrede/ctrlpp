#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_ARM_FACTORY_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_ARM_FACTORY_H

#include "arm.h"
#include "registry.h"
#include "hpipm_arm.h"
#include "json_writer.h"

#include <memory>

namespace ctrlpp::bench::qp_solvers
{

auto make_arm(const arm_spec& arm, tier level, hpipm_dual_signs hpipm_signs) -> std::unique_ptr<qp_arm>;

void describe_arm(json_writer& out, const arm_spec& arm, tier level);

}

#endif
