#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_COMPETITOR_ARMS_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_COMPETITOR_ARMS_H

#include "arm.h"
#include "registry.h"
#include "json_writer.h"

#include <memory>

namespace ctrlpp::bench::qp_solvers
{

/// The leading box rows are taken as variable bounds and the remaining rows as constraint rows.
auto make_qpoases_arm(tier level) -> std::unique_ptr<qp_arm>;

/// The leading box rows are taken as simple bounds; rows with coinciding bounds as equalities.
auto make_daqp_arm(tier level) -> std::unique_ptr<qp_arm>;

/// Rows with coinciding bounds go to the equality block A x = b, every other row to l <= C x <= u.
auto make_proxqp_dense_arm(tier level) -> std::unique_ptr<qp_arm>;

auto make_proxqp_sparse_arm(tier level) -> std::unique_ptr<qp_arm>;

void describe_qpoases(json_writer& out, tier level);

void describe_daqp(json_writer& out, tier level);

void describe_proxqp(json_writer& out, tier level, bool sparse);

}

#endif
