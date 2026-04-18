// Competitive benchmark: ctrlpp::care vs ct::optcon::CARE
// Problem: continuous-time Riccati solve for damped chain-of-integrators systems, size-swept NX in {2, 4, 8}
//
// Warning: ct_optcon has catkin (ROS) heritage. Install via AUR: yay -S control-toolbox

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/control/care.h"

#include <cassert>  // must precede ct_optcon includes; DynamicRiccatiEquation.hpp uses assert() without <cassert>
#include <ct/optcon/lqr/riccati/CARE.hpp>
#include <ct/optcon/lqr/riccati/CARE-impl.hpp>

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
auto build_damped_chain()
{
    // Continuous-time damped chain: A has -0.5 on the diagonal and 1.0 on the superdiagonal.
    // Spectrum is Re(lambda) < 0 so CARE is well-defined.
    Eigen::Matrix<double, int(NX), int(NX)> A = Eigen::Matrix<double, int(NX), int(NX)>::Zero();
    for(std::size_t i = 0; i < NX; ++i)
        A(int(i), int(i)) = -0.5;
    for(std::size_t i = 0; i + 1 < NX; ++i)
        A(int(i), int(i + 1)) = 1.0;

    Eigen::Matrix<double, int(NX), int(NU)> B = Eigen::Matrix<double, int(NX), int(NU)>::Zero();
    const std::size_t group = NX / NU;
    for(std::size_t j = 0; j < NU; ++j)
    {
        const std::size_t last_row = std::min((j + 1) * group, NX) - 1;
        B(int(last_row), int(j)) = 1.0;
    }

    Eigen::Matrix<double, int(NX), int(NX)> Q = Eigen::Matrix<double, int(NX), int(NX)>::Identity();
    Eigen::Matrix<double, int(NU), int(NU)> R = 0.1 * Eigen::Matrix<double, int(NU), int(NU)>::Identity();

    return std::tuple{A, B, Q, R};
}

template <std::size_t NX, std::size_t NU>
void run_size_sweep(ankerl::nanobench::Bench& bench, const char* label_ctrlpp, const char* label_ct)
{
    auto [A, B, Q, R] = build_damped_chain<NX, NU>();

    // ctrlpp warm-up
    auto P_ctrlpp = ctrlpp::care<double, NX, NU>(A, B, Q, R);

    // ct::optcon::CARE setup
    ct::optcon::CARE<NX, NU> ct_care;
    typename ct::optcon::CARE<NX, NU>::state_matrix_t A_ct = A;
    typename ct::optcon::CARE<NX, NU>::control_gain_matrix_t B_ct = B;
    typename ct::optcon::CARE<NX, NU>::state_matrix_t Q_ct = Q;
    typename ct::optcon::CARE<NX, NU>::control_matrix_t R_ct = R;

    // Warm up
    ct_care.computeSteadyStateRiccatiMatrix(Q_ct, R_ct, A_ct, B_ct);

    bench.run(label_ctrlpp,
              [&]
              {
                  auto P = ctrlpp::care<double, NX, NU>(A, B, Q, R);
                  ankerl::nanobench::doNotOptimizeAway(P);
              })
        .run(label_ct,
             [&]
             {
                 auto P_ct = ct_care.computeSteadyStateRiccatiMatrix(Q_ct, R_ct, A_ct, B_ct);
                 ankerl::nanobench::doNotOptimizeAway(P_ct);
             });
}

} // namespace

int main()
{
    ankerl::nanobench::Bench bench;
    bench.title("CARE: ctrlpp vs ct_optcon (size sweep)")
        .warmup(50)
        .minEpochIterations(100)
        .performanceCounters(true)
        .relative(true);

    run_size_sweep<2, 1>(bench, "ctrlpp::care NX=2", "ct::optcon::CARE NX=2");
    run_size_sweep<4, 2>(bench, "ctrlpp::care NX=4", "ct::optcon::CARE NX=4");
    run_size_sweep<8, 2>(bench, "ctrlpp::care NX=8", "ct::optcon::CARE NX=8");

    std::ofstream csv("bench_care_vs_ct.csv");
    bench.render(comma_csv_tpl, csv);
}
