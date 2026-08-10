// Competitive benchmark: ctrlpp::dare vs ct::optcon::DARE
// Problem: discrete-time Riccati solve for chain-of-integrators systems, size-swept
// NX in {2, 4, 6, 8, 12}. NU scales with NX.
//
// Warning: ct_optcon has catkin (ROS) heritage. Install via AUR: yay -S control-toolbox

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "riccati_problem.h"

#include "bench_csv.h"
#include "bench_construct.h"

#include "ctrlpp/control/dare.h"

#include <cassert>  // must precede ct_optcon includes; DynamicRiccatiEquation.hpp uses assert() without <cassert>
#include <ct/optcon/lqr/riccati/DARE.hpp>
#include <ct/optcon/lqr/riccati/DARE-impl.hpp>

#include <Eigen/Dense>

#include <cstddef>
#include <fstream>

namespace
{

using ctrlpp::bench::build_integrator_chain;
using ctrlpp::bench::dare_relative_residual;
using ctrlpp::bench::riccati_plant;

constexpr double sample_period = 0.05;

// ct's DARE writes a gain into a caller-owned matrix alongside the solution it
// returns, so the arm owns that matrix rather than declaring one per call.
template <std::size_t NX, std::size_t NU>
class ct_dare_arm
{
public:
    explicit ct_dare_arm(const riccati_plant<NX, NU>& plant)
        : m_solver{}, m_A{plant.A}, m_Q{plant.Q}, m_R{plant.R}, m_K{}, m_B{plant.B}
    {
    }

    Eigen::Matrix<double, int(NX), int(NX)> solve()
    {
        return m_solver.computeSteadyStateRiccatiMatrix(m_Q, m_R, m_A, m_B, m_K);
    }

private:
    ct::optcon::DARE<NX, NU>                                 m_solver;
    typename ct::optcon::DARE<NX, NU>::state_matrix_t        m_A;
    typename ct::optcon::DARE<NX, NU>::state_matrix_t        m_Q;
    typename ct::optcon::DARE<NX, NU>::control_matrix_t      m_R;
    typename ct::optcon::DARE<NX, NU>::control_feedback_t    m_K;
    typename ct::optcon::DARE<NX, NU>::control_gain_matrix_t m_B;
};

constexpr char const* deviation_metric = "max abs entrywise deviation of the two arms' Riccati solutions P";
constexpr char const* residual_metric = "relative residual of this arm's own discrete Riccati solution";

// The own-criterion rows re-run the identical arm, so their timing columns are a
// second sample of the same work rather than a different workload. Each residual
// is taken on that solver's own solution: routing one solver's P through the
// other's acceptance expression would judge it by a criterion it never met.
template <std::size_t NX, std::size_t NU>
void emit_rows(ankerl::nanobench::Bench& bench, const riccati_plant<NX, NU>& plant, ct_dare_arm<NX, NU>& ct_arm,
               const char* label_ctrlpp, const char* label_ct,
               const Eigen::Matrix<double, int(NX), int(NX)>& P_ctrlpp,
               const Eigen::Matrix<double, int(NX), int(NX)>& P_ct)
{
    auto solve_ctrlpp = [&]
    {
        auto P = ctrlpp::dare<double, NX, NU>(plant.A, plant.B, plant.Q, plant.R);
        ankerl::nanobench::doNotOptimizeAway(P);
    };
    auto solve_ct = [&]
    {
        auto P = ct_arm.solve();
        ankerl::nanobench::doNotOptimizeAway(P);
    };

    ctrlpp::bench::report_accuracy(bench, deviation_metric, (P_ctrlpp - P_ct).cwiseAbs().maxCoeff());
    bench.run(label_ctrlpp, solve_ctrlpp).run(label_ct, solve_ct);
    ctrlpp::bench::run_own_criterion_pair(bench, residual_metric, label_ctrlpp,
                                          dare_relative_residual<NX, NU>(plant, P_ctrlpp), solve_ctrlpp, label_ct,
                                          dare_relative_residual<NX, NU>(plant, P_ct), solve_ct);
}

template <std::size_t NX, std::size_t NU>
void run_size_sweep(ankerl::nanobench::Bench& bench, const char* label_ctrlpp, const char* label_ct)
{
    const riccati_plant<NX, NU> plant = build_integrator_chain<NX, NU>(sample_period);
    const Eigen::Matrix<double, int(NX), int(NX)> P_ctrlpp =
        ctrlpp::bench::built_or_exit(ctrlpp::dare<double, NX, NU>(plant.A, plant.B, plant.Q, plant.R), label_ctrlpp).P;
    ct_dare_arm<NX, NU> ct_arm{plant};
    const Eigen::Matrix<double, int(NX), int(NX)> P_ct = ct_arm.solve();
    emit_rows<NX, NU>(bench, plant, ct_arm, label_ctrlpp, label_ct, P_ctrlpp, P_ct);
}

}

int main(int argc, char** argv)
{
    ankerl::nanobench::Bench bench;
    bench.title("DARE: ctrlpp vs ct_optcon (size sweep)")
        .warmup(50)
        .minEpochIterations(100)
        .performanceCounters(true)
        .relative(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    // The sweep stops at NX = 15. Above it ctrlpp's acceptance check holds an
    // NX(NX+1)/2 square operator as a fixed-size Eigen object, which is 147,968
    // bytes at NX = 16 and 1,729,800 bytes at NX = 30 against Eigen's 131,072
    // byte stack-allocation limit, so the higher rungs cannot be instantiated
    // against the library as shipped.
    run_size_sweep<2, 1>(bench,  "ctrlpp::dare NX=2",  "ct::optcon::DARE NX=2");
    run_size_sweep<4, 2>(bench,  "ctrlpp::dare NX=4",  "ct::optcon::DARE NX=4");
    run_size_sweep<6, 2>(bench,  "ctrlpp::dare NX=6",  "ct::optcon::DARE NX=6");
    run_size_sweep<8, 2>(bench,  "ctrlpp::dare NX=8",  "ct::optcon::DARE NX=8");
    run_size_sweep<12, 3>(bench, "ctrlpp::dare NX=12", "ct::optcon::DARE NX=12");

    std::ofstream csv("bench_dare_vs_ct.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
