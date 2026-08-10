// Competitive benchmark: ctrlpp::lqr_gain (discrete, through the DARE) vs
// ct::optcon::LQR (continuous, through the CARE), on one sampled
// chain-of-integrators plant. Size-swept NX in {2, 4, 6, 8, 12}. NU scales with NX.
//
// The two arms solve different equations on the same matrices, so the gain
// deviation published here is not a rounding figure: it decides whether the
// pairing is a comparison at all. The certificate asks the single question both
// gains must answer -- whether the gain stabilizes the sampled plant it would be
// applied to -- and a radius at or above one says that arm's gain does not.
//
// Warning: ct_optcon has catkin (ROS) heritage. Install via AUR: yay -S control-toolbox

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ct_lqr_arm.h"
#include "riccati_problem.h"

#include "bench_csv.h"
#include "bench_construct.h"

#include "ctrlpp/control/lqr.h"

#include <Eigen/Dense>

#include <cstddef>
#include <fstream>

namespace
{

using ctrlpp::bench::build_integrator_chain;
using ctrlpp::bench::closed_loop_spectral_radius;
using ctrlpp::bench::ct_lqr_arm;
using ctrlpp::bench::riccati_plant;

constexpr double sample_period = 0.05;

constexpr char const* deviation_metric = "max abs entrywise deviation of the two arms' gains K";
constexpr char const* radius_metric =
    "spectral radius of the discrete closed loop this arm's gain forms (below one certifies stability)";

// No own-criterion error row exists here: the arms satisfy different equations,
// so there is no single residual both are entitled to be judged by, and routing
// one arm's gain through the other's equation would score it against a criterion
// it never targeted. The radius is admissible for both because it is a property
// of the sampled plant they are both handed, not of either solver's equation.
template <std::size_t NX, std::size_t NU>
void emit_rows(ankerl::nanobench::Bench& bench, const riccati_plant<NX, NU>& plant, ct_lqr_arm<NX, NU>& ct_arm,
               const char* label_ctrlpp, const char* label_ct,
               const Eigen::Matrix<double, int(NU), int(NX)>& K_ctrlpp,
               const Eigen::Matrix<double, int(NU), int(NX)>& K_ct)
{
    auto solve_ctrlpp = [&]
    {
        auto K = ctrlpp::lqr_gain<double, NX, NU>(plant.A, plant.B, plant.Q, plant.R);
        ankerl::nanobench::doNotOptimizeAway(K);
    };
    auto solve_ct = [&]
    {
        ct_arm.solve();
        ankerl::nanobench::doNotOptimizeAway(ct_arm.gain());
    };

    ctrlpp::bench::report_accuracy(bench, deviation_metric, (K_ctrlpp - K_ct).cwiseAbs().maxCoeff());
    bench.run(label_ctrlpp, solve_ctrlpp).run(label_ct, solve_ct);
    ctrlpp::bench::run_certificate_pair(bench, radius_metric, label_ctrlpp,
                                        closed_loop_spectral_radius<NX, NU>(plant, K_ctrlpp), solve_ctrlpp, label_ct,
                                        closed_loop_spectral_radius<NX, NU>(plant, K_ct), solve_ct);
}

template <std::size_t NX, std::size_t NU>
void run_size_sweep(ankerl::nanobench::Bench& bench, const char* label_ctrlpp, const char* label_ct)
{
    const riccati_plant<NX, NU> plant = build_integrator_chain<NX, NU>(sample_period);
    const Eigen::Matrix<double, int(NU), int(NX)> K_ctrlpp = ctrlpp::bench::built_or_exit(
        ctrlpp::lqr_gain<double, NX, NU>(plant.A, plant.B, plant.Q, plant.R), label_ctrlpp);
    ct_lqr_arm<NX, NU> ct_arm{plant};
    ct_arm.solve();
    const Eigen::Matrix<double, int(NU), int(NX)> K_ct = ct_arm.gain();
    emit_rows<NX, NU>(bench, plant, ct_arm, label_ctrlpp, label_ct, K_ctrlpp, K_ct);
}

}

int main(int argc, char** argv)
{
    ankerl::nanobench::Bench bench;
    bench.title("Discrete LQR: ctrlpp::lqr_gain (DARE) vs ct::optcon::LQR (CARE) (size sweep)")
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
    run_size_sweep<2, 1>(bench,  "ctrlpp::lqr_gain NX=2",  "ct::optcon::LQR NX=2");
    run_size_sweep<4, 2>(bench,  "ctrlpp::lqr_gain NX=4",  "ct::optcon::LQR NX=4");
    run_size_sweep<6, 2>(bench,  "ctrlpp::lqr_gain NX=6",  "ct::optcon::LQR NX=6");
    run_size_sweep<8, 2>(bench,  "ctrlpp::lqr_gain NX=8",  "ct::optcon::LQR NX=8");
    run_size_sweep<12, 3>(bench, "ctrlpp::lqr_gain NX=12", "ct::optcon::LQR NX=12");

    std::ofstream csv("bench_lqr_vs_ct.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
