// Competitive benchmark: ctrlpp::mpc vs TinyMPC (ADMM, embedded LMPC)
// Problem: discrete double-integrator NX=4 NU=2 horizon N=10, box-bounded.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "lmpc/double_integrator.h"

#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/model/state_space.h"

#include <tinympc/tiny_api.hpp>
#include <tinympc/types.hpp>

#include <Eigen/Dense>

#include <cstdint>
#include <fstream>

namespace
{

constexpr char const* comma_csv_tpl = R"TEMPLATE(
"title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

}

int main()
{
    namespace problems = ctrlpp::bench::problems::lmpc;

    constexpr std::size_t NX = 4;
    constexpr std::size_t NU = 2;
    constexpr int N = 10;

    auto sys = problems::make_double_integrator_4_2_state_space();
    auto cfg = problems::make_double_integrator_4_2_config(N);
    auto x0  = problems::double_integrator_4_2_x0_default();

    ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver> ctrlpp_mpc(sys, cfg);
    [[maybe_unused]] auto warm = ctrlpp_mpc.solve(x0);

    // ---- TinyMPC -----------------------------------------------------------
    tinyMatrix A_dyn = sys.A.cast<tinytype>();
    tinyMatrix B_dyn = sys.B.cast<tinytype>();
    tinyMatrix f_dyn = tinyVector::Zero(NX);
    tinyMatrix Q     = cfg.Q.cast<tinytype>();
    tinyMatrix R     = cfg.R.cast<tinytype>();

    TinySolver* tiny = nullptr;
    int rc = tiny_setup(&tiny, A_dyn, B_dyn, f_dyn, Q, R,
                        /*rho=*/static_cast<tinytype>(1.0),
                        static_cast<int>(NX), static_cast<int>(NU), N,
                        /*verbose=*/0);
    if(rc != 0)
    {
        std::fprintf(stderr, "tiny_setup failed: %d\n", rc);
        return 1;
    }

    tinyMatrix x_min_mat = (*cfg.x_min).cast<tinytype>().replicate(1, N + 1);
    tinyMatrix x_max_mat = (*cfg.x_max).cast<tinytype>().replicate(1, N + 1);
    tinyMatrix u_min_mat = (*cfg.u_min).cast<tinytype>().replicate(1, N);
    tinyMatrix u_max_mat = (*cfg.u_max).cast<tinytype>().replicate(1, N);
    tiny_set_bound_constraints(tiny, x_min_mat, x_max_mat, u_min_mat, u_max_mat);

    tinyMatrix x_ref = tinyMatrix::Zero(NX, N + 1);
    tinyMatrix u_ref = tinyMatrix::Zero(NU, N);
    tiny_set_x_ref(tiny, x_ref);
    tiny_set_u_ref(tiny, u_ref);

    tinyVector tiny_x0 = x0.cast<tinytype>();
    tiny_set_x0(tiny, tiny_x0);
    tiny_solve(tiny);

    // ---- Benchmark ---------------------------------------------------------
    ankerl::nanobench::Bench bench;
    bench.title("LMPC: ctrlpp vs TinyMPC")
        .warmup(50)
        .minEpochIterations(100)
        .performanceCounters(true)
        .relative(true)
        .run("ctrlpp::mpc::solve",
             [&]
             {
                 auto u = ctrlpp_mpc.solve(x0);
                 ankerl::nanobench::doNotOptimizeAway(u);
             })
        .run("tiny_solve",
             [&]
             {
                 tiny_set_x0(tiny, tiny_x0);
                 tiny_solve(tiny);
                 ankerl::nanobench::doNotOptimizeAway(tiny->solution);
             });

    std::ofstream csv("bench_lmpc_vs_tinympc.csv");
    bench.render(comma_csv_tpl, csv);
}
