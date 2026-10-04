#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_BINDING_GATE_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_BINDING_GATE_H

#include "hpipm_arm.h"
#include "json_writer.h"

#include <string>
#include <vector>
#include <cstdint>

namespace ctrlpp::bench::qp_solvers
{

/// Over the closed-loop ticks of a cell's tight variant: the median of the binding fraction of its
/// inequality rows, and the largest fraction of ticks on which any single state-box row binds.
struct gate_cell
{
    std::string cell;
    double median_binding_fraction;
    double max_state_row_step_fraction;
    std::int32_t ticks_without_reference;
};

struct gate_family
{
    std::string family;
    std::vector<gate_cell> cells;
    bool passed;
};

auto binding_gate(hpipm_dual_signs signs) -> std::vector<gate_family>;

void write_gate(json_writer& out, const std::vector<gate_family>& families);

}

#endif
