// Competitive benchmark: ctrlpp::care vs ct::optcon::CARE
// Problem: continuous-time Riccati solve for damped chain-of-integrators systems,
// size-swept NX in {2, 4, 6, 8, 12, 16, 20, 24, 30}. NU scales with NX.
//
// Warning: ct_optcon has catkin (ROS) heritage. Install via AUR: yay -S control-toolbox

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_csv.h"
#include "bench_construct.h"

#include "ctrlpp/control/care.h"

#include <cassert>  // must precede ct_optcon includes; DynamicRiccatiEquation.hpp uses assert() without <cassert>
#include <ct/optcon/lqr/riccati/CARE.hpp>
#include <ct/optcon/lqr/riccati/CARE-impl.hpp>

#include <Eigen/Dense>

#include <cstddef>
#include <fstream>

namespace
{

template <std::size_t NX, std::size_t NU>
struct damped_chain
{
    Eigen::Matrix<double, int(NX), int(NX)> A;
    Eigen::Matrix<double, int(NX), int(NU)> B;
    Eigen::Matrix<double, int(NX), int(NX)> Q;
    Eigen::Matrix<double, int(NU), int(NU)> R;
};

template <std::size_t NX, std::size_t NU>
damped_chain<NX, NU> build_damped_chain()
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

    return damped_chain<NX, NU>{A, B, Q, R};
}

template <std::size_t NX, std::size_t NU>
class ct_care_arm
{
public:
    explicit ct_care_arm(const damped_chain<NX, NU>& plant)
        : m_solver{}, m_A{plant.A}, m_Q{plant.Q}, m_R{plant.R}, m_B{plant.B}
    {
    }

    Eigen::Matrix<double, int(NX), int(NX)> solve()
    {
        return m_solver.computeSteadyStateRiccatiMatrix(m_Q, m_R, m_A, m_B);
    }

private:
    ct::optcon::CARE<NX, NU>                                m_solver;
    typename ct::optcon::CARE<NX, NU>::state_matrix_t        m_A;
    typename ct::optcon::CARE<NX, NU>::state_matrix_t        m_Q;
    typename ct::optcon::CARE<NX, NU>::control_matrix_t      m_R;
    typename ct::optcon::CARE<NX, NU>::control_gain_matrix_t m_B;
};

// Relative residual of the continuous algebraic Riccati equation: the residual
// norm over the sum of the norms of the terms that cancel to form it, so the
// figure is dimensionless and comparable across the size sweep.
template <std::size_t NX, std::size_t NU>
double care_relative_residual(const damped_chain<NX, NU>& plant, const Eigen::Matrix<double, int(NX), int(NX)>& P)
{
    const Eigen::Matrix<double, int(NX), int(NX)> cross = plant.A.transpose() * P + P * plant.A;
    const Eigen::Matrix<double, int(NX), int(NX)> quad =
        P * plant.B * plant.R.inverse() * plant.B.transpose() * P;
    return (cross - quad + plant.Q).norm() / (cross.norm() + quad.norm() + plant.Q.norm());
}

constexpr char const* deviation_metric = "max abs entrywise deviation of the two arms' Riccati solutions P";
constexpr char const* residual_metric = "relative residual of this arm's own continuous Riccati solution";

// The own-criterion rows re-run the identical arm, so their timing columns are a
// second sample of the same work rather than a different workload. Each residual
// is taken on that solver's own solution: routing one solver's P through the
// other's acceptance expression would judge it by a criterion it never met.
template <std::size_t NX, std::size_t NU>
void emit_rows(ankerl::nanobench::Bench& bench, const damped_chain<NX, NU>& plant, ct_care_arm<NX, NU>& ct_arm,
               const char* label_ctrlpp, const char* label_ct,
               const Eigen::Matrix<double, int(NX), int(NX)>& P_ctrlpp,
               const Eigen::Matrix<double, int(NX), int(NX)>& P_ct)
{
    auto solve_ctrlpp = [&]
    {
        auto P = ctrlpp::care<double, NX, NU>(plant.A, plant.B, plant.Q, plant.R);
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
                                          care_relative_residual<NX, NU>(plant, P_ctrlpp), solve_ctrlpp, label_ct,
                                          care_relative_residual<NX, NU>(plant, P_ct), solve_ct);
}

template <std::size_t NX, std::size_t NU>
void run_size_sweep(ankerl::nanobench::Bench& bench, const char* label_ctrlpp, const char* label_ct)
{
    const damped_chain<NX, NU> plant = build_damped_chain<NX, NU>();
    const Eigen::Matrix<double, int(NX), int(NX)> P_ctrlpp =
        ctrlpp::bench::built_or_exit(ctrlpp::care<double, NX, NU>(plant.A, plant.B, plant.Q, plant.R), label_ctrlpp)
            .P;
    ct_care_arm<NX, NU> ct_arm{plant};
    const Eigen::Matrix<double, int(NX), int(NX)> P_ct = ct_arm.solve();
    emit_rows<NX, NU>(bench, plant, ct_arm, label_ctrlpp, label_ct, P_ctrlpp, P_ct);
}

} // namespace

int main(int argc, char** argv)
{
    ankerl::nanobench::Bench bench;
    bench.title("CARE: ctrlpp vs ct_optcon (size sweep)")
        .warmup(50)
        .minEpochIterations(100)
        .performanceCounters(true)
        .relative(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    run_size_sweep<2, 1>(bench,  "ctrlpp::care NX=2",  "ct::optcon::CARE NX=2");
    run_size_sweep<4, 2>(bench,  "ctrlpp::care NX=4",  "ct::optcon::CARE NX=4");
    run_size_sweep<6, 2>(bench,  "ctrlpp::care NX=6",  "ct::optcon::CARE NX=6");
    run_size_sweep<8, 2>(bench,  "ctrlpp::care NX=8",  "ct::optcon::CARE NX=8");
    run_size_sweep<12, 3>(bench, "ctrlpp::care NX=12", "ct::optcon::CARE NX=12");
    run_size_sweep<16, 4>(bench, "ctrlpp::care NX=16", "ct::optcon::CARE NX=16");
    run_size_sweep<20, 5>(bench, "ctrlpp::care NX=20", "ct::optcon::CARE NX=20");
    run_size_sweep<24, 6>(bench, "ctrlpp::care NX=24", "ct::optcon::CARE NX=24");
    run_size_sweep<30, 6>(bench, "ctrlpp::care NX=30", "ct::optcon::CARE NX=30");

    std::ofstream csv("bench_care_vs_ct.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
