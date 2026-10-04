#include "probe.h"
#include "driver.h"
#include "describe.h"
#include "blocks.h"
#include "binding_gate.h"
#include "driver_options.h"

#include <cstdio>
#include <format>
#include <fstream>
#include <algorithm>

namespace ctrlpp::bench::qp_solvers
{

namespace
{

auto output_path(const driver_options& options, std::string_view suffix) -> std::string
{
    return std::format("{}/{}_rep{}{}", options.out_dir, name_of(options.group), options.rep, suffix);
}

auto run_gate(const driver_options& options, const hpipm_settlement& hpipm) -> int
{
    const std::vector<gate_family> families = binding_gate(hpipm.best);
    std::ofstream file(options.gate_path);
    json_writer out(file);
    write_gate(out, families);
    for(const gate_family& family : families)
        for(const gate_cell& cell : family.cells)
            std::printf("binding gate %-8s %-13s G1 %.4f G2 %.4f%s %s\n", family.family.c_str(), cell.cell.c_str(),
                        cell.median_binding_fraction, cell.max_state_row_step_fraction,
                        cell.ticks_without_reference > 0 ? std::format(" ({} ticks without a reference)", cell.ticks_without_reference).c_str() : "",
                        family.passed ? "family passed" : "family refused");
    return file.good() ? 0 : 1;
}

auto run_probe_stage(const driver_options& options, const hpipm_settlement& hpipm) -> bool
{
    const std::vector<probe_result> results = run_probes(options.group, hpipm);
    std::ofstream file(output_path(options, "_probes.json"));
    json_writer out(file);
    write_probes(out, results, hpipm);
    bool all = true;
    for(const probe_result& result : results)
    {
        std::fprintf(stderr, "probe %-28s %s max_dual_error %.3e\n", result.arm.c_str(), result.passed ? "passed" : "FAILED", result.max_dual_error);
        all = all && result.passed;
    }
    return all;
}

}

auto run_driver(cohort group, int argc, char** argv) -> int
{
    const auto options = parse_options(group, argc, argv);
    if(!options)
        return 2;
    if(!options->describe_path.empty())
        return write_description(options->describe_path) ? 0 : 1;
    const hpipm_settlement hpipm = settle_hpipm_signs();
    if(!options->gate_path.empty())
        return run_gate(*options, hpipm);
    if(options->run_id.empty())
    {
        std::fprintf(stderr, "qp_solvers_%s: a measured run needs --run-id\n", name_of(group).data());
        return 2;
    }
    if(!run_probe_stage(*options, hpipm) && options->on_probe_failure == "abort")
    {
        std::fprintf(stderr, "qp_solvers_%s: a convention probe failed; nothing was measured\n", name_of(group).data());
        return 3;
    }
    record_writer records(output_path(*options, ".csv"), options->run_id);
    run_blocks(*options, hpipm.best, records);
    return records.good() ? 0 : 1;
}

}
