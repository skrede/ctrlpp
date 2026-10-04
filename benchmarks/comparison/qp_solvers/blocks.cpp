#include "cases.h"
#include "blocks.h"
#include "sequence.h"
#include "frequency_log.h"

#include "qp/mpc_tracking_cells.h"

#include <cstdio>
#include <format>
#include <memory>
#include <optional>
#include <algorithm>

namespace ctrlpp::bench::qp_solvers
{

namespace
{

namespace qp = problems::qp;

struct case_forms
{
    std::unique_ptr<problem_case> dense;
    std::unique_ptr<problem_case> sparse;
    const problem_case& source;
};

auto expressible(const arm_spec& arm, std::int32_t horizon) -> bool
{
    return arm.kind != arm_kind::tinympc || horizon >= qp::long_horizon_min_steps;
}

auto block_arms(cohort group, tier level, std::int32_t horizon, std::int32_t rep) -> std::vector<arm_spec>
{
    std::vector<arm_spec> arms;
    for(const arm_spec& arm : arm_registry)
        if(arm.group == group && runs_at(arm, level) && expressible(arm, horizon))
            arms.push_back(arm);
    if(!arms.empty())
        std::ranges::rotate(arms, arms.begin() + rep % static_cast<std::int32_t>(arms.size()));
    return arms;
}

void run_tiers(const driver_options& options, const sequence_block& base, const sequence_tools& tools)
{
    for(const tier level : tiers)
    {
        if(!admits(options.tiers, name_of(level)))
            continue;
        sequence_block block = base;
        block.level = level;
        std::fprintf(stderr, "%s rep %d  %s %s %s %s\n", name_of(block.group).data(), block.rep, block.cell.c_str(),
                     block.variant.data(), name_of(block.kind).data(), name_of(level).data());
        for(const arm_spec& arm : block_arms(block.group, level, block.N, block.rep))
            if(admits(options.arms, arm.id))
                run_sequence(block, arm, tools);
    }
}

void run_regimes(const driver_options& options, const case_forms& forms, const sequence_block& base, hpipm_dual_signs signs,
                 const sequence_tools& tools)
{
    const bool shaped = base.nx == 0;
    for(const regime kind : regimes)
    {
        if(!admits(options.regimes, name_of(kind)) || (shaped && kind != regime::fixed))
            continue;
        const std::int32_t seed = kind == regime::churning ? options.rep : -1;
        const auto references = reference_sequence({*forms.dense, forms.sparse.get()}, kind, seed, options.steps, signs);
        sequence_block block{base.group, base.source, base.cell, base.variant, base.nx, base.nu, base.N, kind, base.level,
                             base.rep, base.steps, references};
        run_tiers(options, block, tools);
    }
}

auto forms_of(cohort group, const std::optional<qp::mpc_sweep_cell>& cell, const qp::mpc_variant& variant) -> case_forms
{
    if(!cell)
    {
        auto shaped = make_shaped_case();
        const problem_case& source = *shaped;
        return {std::move(shaped), nullptr, source};
    }
    auto dense = make_condensed_case(*cell, variant);
    auto sparse = make_sparse_case(*cell, variant);
    const problem_case& source = group == cohort::dense ? *dense : *sparse;
    return {std::move(dense), std::move(sparse), source};
}

void run_case(const driver_options& options, const std::optional<qp::mpc_sweep_cell>& cell, const qp::mpc_variant& variant,
              hpipm_dual_signs signs, const sequence_tools& tools)
{
    tools.frequencies.sample(cell ? qp::cell_id(*cell) : std::string(shaped_cell_id));
    const case_forms forms = forms_of(options.group, cell, variant);
    const std::vector<reference_step> none;
    const qp::mpc_sweep_cell dims = cell.value_or(qp::mpc_sweep_cell{0, 0, 0});
    const std::string id = cell ? qp::cell_id(*cell) : std::string(shaped_cell_id);
    const sequence_block base{options.group, forms.source, id, variant.name, dims.nx, dims.nu, dims.horizon, regime::fixed,
                              tier::mpc, options.rep, options.steps, none};
    run_regimes(options, forms, base, signs, tools);
}

auto refused(const driver_options& options, const qp::mpc_sweep_cell& cell, const qp::mpc_variant& variant) -> bool
{
    return variant.state_limit && std::ranges::find(options.refused_families, qp::family_id(cell)) != options.refused_families.end();
}

}

void run_blocks(const driver_options& options, hpipm_dual_signs signs, record_writer& records)
{
    const instruction_counter counter;
    frequency_log frequencies(std::format("{}/{}_rep{}_frequency.txt", options.out_dir, name_of(options.group), options.rep));
    const sequence_tools tools{counter, signs, records, frequencies};
    for(const qp::mpc_sweep_cell& cell : qp::tracking_cells())
        for(const qp::mpc_variant& variant : qp::mpc_variants())
            if(admits(options.cells, qp::cell_id(cell)) && admits(options.variants, variant.name) && !refused(options, cell, variant))
                run_case(options, cell, variant, signs, tools);
    const bool shaped = options.group == cohort::dense && admits(options.cells, shaped_cell_id);
    if(shaped && admits(options.variants, qp::nominal_variant().name))
        run_case(options, std::nullopt, qp::nominal_variant(), signs, tools);
}

}
