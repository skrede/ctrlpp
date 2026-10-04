#include "cases.h"
#include "reference.h"
#include "binding_gate.h"

#include "qp/mpc_tracking_cells.h"

#include <algorithm>

namespace ctrlpp::bench::qp_solvers
{

namespace
{

namespace qp = problems::qp;

/// The mean of the two middle values for an even count.
auto median(std::vector<double> values) -> double
{
    if(values.empty())
        return 0.0;
    std::ranges::sort(values);
    const std::size_t half = values.size() / 2;
    return values.size() % 2 == 1 ? values[half] : 0.5 * (values[half - 1] + values[half]);
}

auto state_row_step_fraction(const std::vector<reference_step>& steps, std::size_t first_state_row) -> double
{
    double largest = 0.0;
    std::size_t found = 0;
    std::vector<std::size_t> binding;
    for(const reference_step& step : steps)
    {
        if(!step.found)
            continue;
        ++found;
        binding.resize(step.rows.size(), 0);
        for(std::size_t i = first_state_row; i < step.rows.size(); ++i)
            binding[i] += step.rows[i] != qp::row_state::inactive ? 1 : 0;
    }
    for(std::size_t count : binding)
        largest = std::max(largest, found > 0 ? static_cast<double>(count) / static_cast<double>(found) : 0.0);
    return largest;
}

auto gate_of(const qp::mpc_sweep_cell& cell, hpipm_dual_signs signs) -> gate_cell
{
    const auto dense = make_condensed_case(cell, qp::tight_variant());
    const auto sparse = make_sparse_case(cell, qp::tight_variant());
    const auto steps = reference_sequence({*dense, sparse.get()}, regime::closed_loop, -1, steps_per_sequence, signs);
    std::vector<double> fractions;
    for(const reference_step& step : steps)
        if(step.found)
            fractions.push_back(binding_fraction(step));
    const auto missing = static_cast<std::int32_t>(steps.size() - fractions.size());
    const auto first_state_row = static_cast<std::size_t>(cell.horizon * cell.nu);
    return {qp::cell_id(cell), median(fractions), state_row_step_fraction(steps, first_state_row), missing};
}

auto family_passes(const std::vector<gate_cell>& cells) -> bool
{
    return std::ranges::all_of(cells, [](const gate_cell& c) {
        return c.median_binding_fraction >= binding_median_fraction_min
            && c.max_state_row_step_fraction >= binding_state_row_step_fraction_min;
    });
}

}

auto binding_gate(hpipm_dual_signs signs) -> std::vector<gate_family>
{
    std::vector<gate_family> families;
    for(const qp::mpc_sweep_cell& cell : qp::tracking_cells())
    {
        const std::string family = qp::family_id(cell);
        auto match = std::ranges::find(families, family, &gate_family::family);
        if(match == families.end())
            match = families.insert(families.end(), gate_family{family, {}, false});
        match->cells.push_back(gate_of(cell, signs));
    }
    for(gate_family& family : families)
        family.passed = family_passes(family.cells);
    return families;
}

void write_gate(json_writer& out, const std::vector<gate_family>& families)
{
    out.open();
    for(const gate_family& family : families)
    {
        out.open(family.family);
        out.open("cells");
        for(const gate_cell& cell : family.cells)
        {
            out.open(cell.cell);
            out.number("median_binding_fraction", cell.median_binding_fraction);
            out.number("max_state_row_step_fraction", cell.max_state_row_step_fraction);
            out.close();
        }
        out.close();
        out.flag("family_passed", family.passed);
        out.close();
    }
    out.close();
}

}
