#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_SEQUENCE_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_SEQUENCE_H

#include "cases.h"
#include "registry.h"
#include "hpipm_arm.h"
#include "reference.h"
#include "frequency_log.h"
#include "record_writer.h"
#include "instruction_counter.h"

#include <string>
#include <vector>
#include <cstdint>

namespace ctrlpp::bench::qp_solvers
{

struct sequence_block
{
    cohort group;
    const problem_case& source;
    std::string cell;
    std::string_view variant;
    std::int32_t nx;
    std::int32_t nu;
    std::int32_t N;
    regime kind;
    tier level;
    std::int32_t rep;
    std::int32_t steps;
    const std::vector<reference_step>& references;
};

struct sequence_tools
{
    const instruction_counter& counter;
    hpipm_dual_signs hpipm_signs;
    record_writer& records;
    frequency_log& frequencies;
};

/// A fresh arm, set up on the step-0 tick, then (in the fixed regime, after untimed warm resolves)
/// one timed update and solve per step; every answer is judged after its tick.
void run_sequence(const sequence_block& block, const arm_spec& arm, const sequence_tools& tools);

}

#endif
