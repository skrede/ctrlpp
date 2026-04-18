// Competitive benchmark: ctrlpp::dare vs ct::optcon::DARE
// Problem: discrete-time Riccati solve for chain-of-integrators systems, size-swept NX in {2, 4, 8}
//
// Warning: ct_optcon has catkin (ROS) heritage. Install via AUR: yay -S control-toolbox

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/control/dare.h"

#include <cassert>  // must precede ct_optcon includes; DynamicRiccatiEquation.hpp uses assert() without <cassert>
#include <ct/optcon/lqr/riccati/DARE.hpp>
#include <ct/optcon/lqr/riccati/DARE-impl.hpp>

#include <Eigen/Dense>

#include <fstream>
#include <cstddef>

namespace
{

constexpr char const* comma_csv_tpl = R"TEMPLATE(
"title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

template <std::size_t NX, std::size_t NU>
auto build_chain_of_integrators(double dt)
{
    Eigen::Matrix<double, int(NX), int(NX)> A = Eigen::Matrix<double, int(NX), int(NX)>::Identity();
    for(std::size_t i = 0; i + 1 < NX; ++i)
        A(int(i), int(i + 1)) = dt;

    Eigen::Matrix<double, int(NX), int(NU)> B = Eigen::Matrix<double, int(NX), int(NU)>::Zero();
    const std::size_t group = NX / NU;
    for(std::size_t j = 0; j < NU; ++j)
    {
        const std::size_t last_row = std::min((j + 1) * group, NX) - 1;
        B(int(last_row), int(j)) = dt;
    }

    Eigen::Matrix<double, int(NX), int(NX)> Q = Eigen::Matrix<double, int(NX), int(NX)>::Identity();
    Eigen::Matrix<double, int(NU), int(NU)> R = 0.1 * Eigen::Matrix<double, int(NU), int(NU)>::Identity();

    return std::tuple{A, B, Q, R};
}

template <std::size_t NX, std::size_t NU>
void run_size_sweep(ankerl::nanobench::Bench& bench, const char* label_ctrlpp, const char* label_ct)
{
    constexpr double dt = 0.05;
    auto [A, B, Q, R] = build_chain_of_integrators<NX, NU>(dt);

    // ctrlpp warm-up
    auto P_ctrlpp = ctrlpp::dare<double, NX, NU>(A, B, Q, R);

    // ct::optcon::DARE setup
    ct::optcon::DARE<NX, NU> ct_dare;
    typename ct::optcon::DARE<NX, NU>::state_matrix_t A_ct = A;
    typename ct::optcon::DARE<NX, NU>::control_gain_matrix_t B_ct = B;
    typename ct::optcon::DARE<NX, NU>::state_matrix_t Q_ct = Q;
    typename ct::optcon::DARE<NX, NU>::control_matrix_t R_ct = R;
    typename ct::optcon::DARE<NX, NU>::control_feedback_t K_ct;

    // Warm up
    ct_dare.computeSteadyStateRiccatiMatrix(Q_ct, R_ct, A_ct, B_ct, K_ct);

    bench.run(label_ctrlpp,
              [&]
              {
                  auto P = ctrlpp::dare<double, NX, NU>(A, B, Q, R);
                  ankerl::nanobench::doNotOptimizeAway(P);
              })
        .run(label_ct,
             [&]
             {
                 auto P_ct = ct_dare.computeSteadyStateRiccatiMatrix(Q_ct, R_ct, A_ct, B_ct, K_ct);
                 ankerl::nanobench::doNotOptimizeAway(P_ct);
             });
}

} // namespace

int main()
{
    ankerl::nanobench::Bench bench;
    bench.title("DARE: ctrlpp vs ct_optcon (size sweep)")
        .warmup(50)
        .minEpochIterations(100)
        .performanceCounters(true)
        .relative(true);

    run_size_sweep<2, 1>(bench, "ctrlpp::dare NX=2", "ct::optcon::DARE NX=2");
    run_size_sweep<4, 2>(bench, "ctrlpp::dare NX=4", "ct::optcon::DARE NX=4");
    run_size_sweep<8, 2>(bench, "ctrlpp::dare NX=8", "ct::optcon::DARE NX=8");

    std::ofstream csv("bench_dare_vs_ct.csv");
    bench.render(comma_csv_tpl, csv);
}
