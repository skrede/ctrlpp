#include "record_writer.h"

#include <cmath>
#include <format>

namespace ctrlpp::bench::qp_solvers
{

namespace
{

constexpr std::string_view header =
        "run_id,cohort,cell,variant,nx,nu,N,n,m,regime,tier,arm,polish,rep,seed,step,setup_ns,update_ns,solve_ns,tick_ns,"
        "condense_ns,instr_setup,instr_update,instr_solve,instr_tick,iterations,kkt_stat,kkt_prim,kkt_comp,kkt_max,"
        "native_status,pass,fail_class,allocs_tick,dual_kind,fallback,binding_frac,churn_frac";

auto real(double value) -> std::string
{
    return std::isnan(value) ? std::string("nan") : std::format("{:.17g}", value);
}

}

record_writer::record_writer(const std::string& path, std::string_view run_id)
    : out_(path)
    , run_id_(run_id)
{
    out_ << header << '\n';
}

void record_writer::write(const tick_record& r)
{
    out_ << std::format("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},", run_id_, r.cohort, r.cell, r.variant, r.nx, r.nu,
                        r.N, r.n, r.m, r.regime, r.tier, r.arm, r.polish, r.rep, r.seed, r.step);
    out_ << std::format("{},{},{},{},{},{},{},{},{},{},", r.setup_ns, r.update_ns, r.solve_ns, r.tick_ns, r.condense_ns,
                        r.instr_setup, r.instr_update, r.instr_solve, r.instr_tick, r.iterations);
    out_ << std::format("{},{},{},{},{},{},{},{},{},{},{},{}\n", real(r.kkt_stat), real(r.kkt_prim), real(r.kkt_comp),
                        real(r.kkt_max), r.native_status, r.pass, r.fail_class, r.allocs_tick, r.dual_kind, r.fallback,
                        real(r.binding_frac), real(r.churn_frac));
}

}
