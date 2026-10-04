#include "judge.h"
#include "sequence.h"
#include "tick_timer.h"
#include "arm_factory.h"
#include "dual_recovery.h"

namespace ctrlpp::bench::qp_solvers
{

namespace
{

struct setup_sample
{
    std::int64_t ns;
    std::int64_t instructions;
};

auto timed_setup(qp_arm& arm, const qp_instance& instance, const qp_tick& tick, const instruction_counter& counter) -> setup_sample
{
    const std::int64_t i0 = counter.read();
    const tick_clock::time_point t0 = tick_clock::now();
    arm.setup(instance, tick);
    const tick_clock::time_point t1 = tick_clock::now();
    const std::int64_t i1 = counter.read();
    return {elapsed_ns(t0, t1), counted(i0, i1)};
}

/// The cohort's own time to rebuild a tick's vectors from x0; the sparse layout writes x0 into two
/// bound vectors and has nothing to condense, and the dense-only problem has no x0 at all.
auto timed_fill(const sequence_block& block, const Eigen::VectorXd& x0, qp_tick& tick) -> std::int64_t
{
    const tick_clock::time_point t0 = tick_clock::now();
    block.source.fill(x0, tick);
    const tick_clock::time_point t1 = tick_clock::now();
    if(block.group == cohort::sparse)
        return -1;
    return block.nx == 0 ? 0 : elapsed_ns(t0, t1);
}

auto polish_flag(const arm_spec& arm, tier level) -> std::int32_t
{
    const bool has_polish = arm.kind == arm_kind::argmin_dense_admm || arm.kind == arm_kind::argmin_sparse_admm || arm.kind == arm_kind::osqp;
    return has_polish ? (level == tier::high || arm.forced_polish ? 1 : 0) : -1;
}

struct judged
{
    kkt_terms terms;
    bool passed;
};

auto judge_answer(const sequence_block& block, const arm_spec& arm, problems::qp::dense_program& program, const qp_tick& tick,
                  arm_answer& answer) -> judged
{
    program.q = tick.q;
    program.l = tick.l;
    program.u = tick.u;
    if(arm.duals == dual_kind::recovered)
        answer.y = recover_duals(program, block.source.instance().A_sparse, answer.x, threshold(block.level));
    const kkt_terms terms = judge(program, answer.x, answer.y);
    return {terms, passes(terms, block.level)};
}

auto churn_at(const sequence_block& block, std::int32_t step) -> double
{
    if(block.kind != regime::churning || step == 0)
        return -1.0;
    return churn_fraction(block.references[static_cast<std::size_t>(step - 1)], block.references[static_cast<std::size_t>(step)]);
}

void write_row(const sequence_block& block, const arm_spec& arm, std::int32_t step, const tick_sample& sample,
               std::int64_t condense_ns, const setup_sample& setup, const arm_answer& answer, const judged& verdict,
               record_writer& records)
{
    const qp_instance& instance = block.source.instance();
    const bool first = step == 0;
    const std::int32_t seed = block.kind == regime::churning ? block.rep : -1;
    records.write({name_of(block.group), block.cell, block.variant, block.nx, block.nu, block.N, instance.P.rows(), instance.A.rows(),
                   name_of(block.kind), name_of(block.level), arm.id, polish_flag(arm, block.level), block.rep, seed, step,
                   first ? setup.ns : -1, sample.update_ns, sample.solve_ns, sample.tick_ns, condense_ns,
                   first ? setup.instructions : -1, sample.instr_update, sample.instr_solve, sample.instr_tick, answer.iterations,
                   verdict.terms.stationarity, verdict.terms.primal, verdict.terms.complementarity, verdict.terms.worst,
                   answer.status, verdict.passed ? 1 : 0, fail_class(verdict.passed, answer.claim), sample.allocs,
                   arm.duals == dual_kind::native ? "native" : "recovered", answer.fallback ? 1 : 0,
                   binding_fraction(block.references[static_cast<std::size_t>(step)]), churn_at(block, step)});
}

}

void run_sequence(const sequence_block& block, const arm_spec& arm, const sequence_tools& tools)
{
    const qp_instance& instance = block.source.instance();
    const std::int32_t seed = block.kind == regime::churning ? block.rep : -1;
    problems::qp::dense_program program{instance.P, {}, instance.A, {}, {}};
    qp_tick tick;
    block.source.fill(initial_state(block.source, block.kind, seed, 0), tick);
    const auto solver = make_arm(arm, block.level, tools.hpipm_signs);
    const setup_sample setup = timed_setup(*solver, instance, tick, tools.counter);
    for(std::int32_t w = 0; block.kind == regime::fixed && w < fixed_warmup_resolves; ++w)
    {
        solver->update(tick);
        solver->solve();
    }
    arm_answer answer;
    for(std::int32_t step = 0; step < block.steps; ++step)
    {
        const std::int64_t condense_ns = timed_fill(block, initial_state(block.source, block.kind, seed, step), tick);
        const tick_sample sample = timed_tick(*solver, tick, tools.counter);
        solver->answer(answer);
        const judged verdict = judge_answer(block, arm, program, tick, answer);
        write_row(block, arm, step, sample, condense_ns, setup, answer, verdict, tools.records);
    }
}

}
