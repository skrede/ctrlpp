#include "probe.h"
#include "arm_factory.h"
#include "probe_problems.h"

#include <array>
#include <cmath>
#include <limits>

namespace ctrlpp::bench::qp_solvers
{

namespace
{

auto dual_error(qp_arm& arm, const probe_problem& probe) -> double
{
    arm.setup(probe.instance, probe.tick);
    arm.update(probe.tick);
    arm.solve();
    arm_answer answer;
    arm.answer(answer);
    if(answer.y.size() != probe.exact.y.size() || !answer.y.allFinite())
        return std::numeric_limits<double>::infinity();
    return (answer.y - probe.exact.y).lpNorm<Eigen::Infinity>();
}

auto probe_tier(const arm_spec& arm) -> tier
{
    return arm.runs_high ? tier::high : tier::mpc;
}

auto probe_of(const arm_spec& arm, const probe_problem& probe, const hpipm_settlement& hpipm) -> probe_result
{
    if(arm.kind == arm_kind::hpipm)
        return {std::string(arm.id), hpipm.signs.has_value(), hpipm.error};
    const auto solver = make_arm(arm, probe_tier(arm), hpipm.best);
    const double error = dual_error(*solver, probe);
    return {std::string(arm.id), error <= probe_tolerance, error};
}

}

auto settle_hpipm_signs() -> hpipm_settlement
{
    const probe_problem probe = sparse_probe();
    constexpr std::array<double, 2> signs{1.0, -1.0};
    hpipm_settlement out{std::nullopt, {1.0, 1.0}, std::numeric_limits<double>::infinity(), 0};
    for(const double initial_state : signs)
        for(const double dynamics : signs)
        {
            const hpipm_dual_signs pair{initial_state, dynamics};
            const auto arm = make_hpipm_arm(tier::high, pair);
            const double error = dual_error(*arm, probe);
            out.matching_pairs += error <= probe_tolerance ? 1 : 0;
            if(error < out.error)
                out = {out.signs, pair, error, out.matching_pairs};
        }
    if(out.matching_pairs == 1)
        out.signs = out.best;
    return out;
}

auto run_probes(cohort group, const hpipm_settlement& hpipm) -> std::vector<probe_result>
{
    const probe_problem probe = group == cohort::dense ? dense_probe() : sparse_probe();
    std::vector<probe_result> results;
    for(const arm_spec& arm : arm_registry)
        if(arm.group == group && arm.duals == dual_kind::native)
            results.push_back(probe_of(arm, probe, hpipm));
    return results;
}

void write_probes(json_writer& out, const std::vector<probe_result>& results, const hpipm_settlement& hpipm)
{
    out.open();
    out.open("probes");
    for(const probe_result& result : results)
    {
        out.open(result.arm);
        out.flag("passed", result.passed);
        out.number("max_dual_error", result.max_dual_error);
        out.close();
    }
    out.close();
    out.open("hpipm_signs");
    out.number("initial_state", hpipm.best.initial_state);
    out.number("dynamics", hpipm.best.dynamics);
    out.integer("matching_pairs", hpipm.matching_pairs);
    out.close();
    out.close();
}

}
