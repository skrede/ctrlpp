#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_ARGMIN_ARMS_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_ARGMIN_ARMS_H

#include "arm.h"
#include "registry.h"
#include "json_writer.h"

#include <memory>

namespace ctrlpp::bench::qp_solvers
{

auto make_argmin_dense_admm_arm(tier level, bool forced_polish) -> std::unique_ptr<qp_arm>;

auto make_argmin_active_set_arm(tier level) -> std::unique_ptr<qp_arm>;

auto make_argmin_sparse_admm_arm(tier level, bool forced_polish) -> std::unique_ptr<qp_arm>;

auto make_argmin_active_set_reference() -> std::unique_ptr<qp_arm>;

void describe_argmin_dense_admm(json_writer& out, tier level, bool forced_polish);

void describe_argmin_sparse_admm(json_writer& out, tier level, bool forced_polish);

void describe_argmin_active_set(json_writer& out, tier level);

}

#endif
