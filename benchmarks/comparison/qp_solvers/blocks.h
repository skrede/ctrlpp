#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_BLOCKS_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_BLOCKS_H

#include "hpipm_arm.h"
#include "record_writer.h"
#include "driver_options.h"

namespace ctrlpp::bench::qp_solvers
{

/// The repetition's blocks in registry order of (cell, variant, regime, tier), skipping the tight variant
/// of a refused family; within a block the arms that run at the tier and can express the cell, in
/// registry order rotated left by the repetition index.
void run_blocks(const driver_options& options, hpipm_dual_signs signs, record_writer& records);

}

#endif
