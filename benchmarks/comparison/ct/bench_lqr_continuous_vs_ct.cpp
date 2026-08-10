// Competitive benchmark: ctrlpp continuous-time LQR (via care + K = R^-1 B^T P) vs ct::optcon::LQR
// Problem: identical to bench_lqr_vs_ct.cpp (integrator dynamics, NX=4, NU=2, dt=0.1)
//
// ct::optcon::LQR::compute is a continuous-time LQR built on ct::optcon::CARE::solve
// (Eigen::RealSchur + LAPACK dtrsen_ reorder). bench_lqr_vs_ct compares against
// ctrlpp::lqr_gain which chains through DARE -- apples to oranges. This benchmark
// chains ctrlpp::care for the continuous-time Riccati and computes K directly, giving
// a straight apples-to-apples race.
//
// Warning: ct_optcon has catkin (ROS) heritage. Install via AUR: yay -S control-toolbox

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "bench_csv.h"
#include "bench_construct.h"

#include "ctrlpp/control/lqr.h"
#include "ctrlpp/control/care.h"

#include <cassert>  // must precede ct_optcon includes; DynamicRiccatiEquation.hpp uses assert() without <cassert>
#include <ct/core/types/StateVector.h>
#include <ct/core/types/ControlVector.h>
#include <ct/optcon/lqr/riccati/CARE.hpp>
#include <ct/optcon/lqr/riccati/CARE-impl.hpp>
#include <ct/optcon/lqr/LQR.hpp>
#include <ct/optcon/lqr/LQR-impl.hpp>

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

// ct's LQR writes its gain into a caller-owned matrix, so the arm holds that
// matrix rather than returning one: a return by value inside the timed region
// would charge the ct arm for a copy the ctrlpp arm never makes.
template <std::size_t NX, std::size_t NU>
class ct_lqr_arm
{
public:
    explicit ct_lqr_arm(const damped_chain<NX, NU>& plant)
        : m_lqr{}, m_K{}, m_B{plant.B}, m_A{plant.A}, m_Q{plant.Q}, m_R{plant.R}
    {
    }

    void solve()
    {
        m_lqr.compute(m_Q, m_R, m_A, m_B, m_K);
    }

    const Eigen::Matrix<double, int(NU), int(NX)>& gain() const
    {
        return m_K;
    }

private:
    ct::optcon::LQR<NX, NU>                            m_lqr;
    Eigen::Matrix<double, int(NU), int(NX)>            m_K;
    Eigen::Matrix<double, int(NX), int(NU)>            m_B;
    typename ct::optcon::LQR<NX, NU>::state_matrix_t   m_A;
    typename ct::optcon::LQR<NX, NU>::state_matrix_t   m_Q;
    typename ct::optcon::LQR<NX, NU>::control_matrix_t m_R;
};

template <std::size_t NX, std::size_t NU>
double closed_loop_abscissa(const damped_chain<NX, NU>& plant, const Eigen::Matrix<double, int(NU), int(NX)>& K)
{
    const Eigen::Matrix<double, int(NX), int(NX)> closed_loop = plant.A - plant.B * K;
    const Eigen::EigenSolver<Eigen::Matrix<double, int(NX), int(NX)>> spectrum(closed_loop, false);
    return spectrum.eigenvalues().real().maxCoeff();
}

constexpr char const* deviation_metric = "max abs entrywise deviation of the two arms' gains K";
constexpr char const* abscissa_metric = "closed-loop spectral abscissa of this arm's own gain";

// A gain deviation alone cannot say the two arms solved the same problem; the
// abscissa is what shows a gain actually stabilizes the plant, and unlike the
// deviation it is defined for one arm alone. The own-criterion rows re-run the
// identical arm, so their timing columns are a second sample of the same work.
template <std::size_t NX, std::size_t NU>
void emit_rows(ankerl::nanobench::Bench& bench, const damped_chain<NX, NU>& plant, ct_lqr_arm<NX, NU>& ct_arm,
               const char* label_ctrlpp, const char* label_ct,
               const Eigen::Matrix<double, int(NU), int(NX)>& K_ctrlpp,
               const Eigen::Matrix<double, int(NU), int(NX)>& K_ct)
{
    auto solve_ctrlpp = [&]
    {
        auto K = ctrlpp::lqr_gain_continuous<double, NX, NU>(plant.A, plant.B, plant.Q, plant.R);
        ankerl::nanobench::doNotOptimizeAway(K);
    };
    auto solve_ct = [&]
    {
        ct_arm.solve();
        ankerl::nanobench::doNotOptimizeAway(ct_arm.gain());
    };

    ctrlpp::bench::report_accuracy(bench, deviation_metric, (K_ctrlpp - K_ct).cwiseAbs().maxCoeff());
    bench.run(label_ctrlpp, solve_ctrlpp).run(label_ct, solve_ct);
    ctrlpp::bench::run_own_criterion_pair(bench, abscissa_metric, label_ctrlpp,
                                          closed_loop_abscissa<NX, NU>(plant, K_ctrlpp), solve_ctrlpp, label_ct,
                                          closed_loop_abscissa<NX, NU>(plant, K_ct), solve_ct);
}

template <std::size_t NX, std::size_t NU>
void run_size_sweep(ankerl::nanobench::Bench& bench, const char* label_ctrlpp, const char* label_ct)
{
    const damped_chain<NX, NU> plant = build_damped_chain<NX, NU>();
    const Eigen::Matrix<double, int(NU), int(NX)> K_ctrlpp = ctrlpp::bench::built_or_exit(
        ctrlpp::lqr_gain_continuous<double, NX, NU>(plant.A, plant.B, plant.Q, plant.R), label_ctrlpp);
    ct_lqr_arm<NX, NU> ct_arm{plant};
    ct_arm.solve();
    const Eigen::Matrix<double, int(NU), int(NX)> K_ct = ct_arm.gain();
    emit_rows<NX, NU>(bench, plant, ct_arm, label_ctrlpp, label_ct, K_ctrlpp, K_ct);
}

} // namespace

int main(int argc, char** argv)
{
    ankerl::nanobench::Bench bench;
    bench.title("Continuous LQR: ctrlpp::lqr_gain_continuous vs ct::optcon::LQR (size sweep)")
        .warmup(50)
        .minEpochIterations(100)
        .performanceCounters(true)
        .relative(true);
    ctrlpp::bench::apply_smoke_switch(bench, argc, argv);

    run_size_sweep<2, 1>(bench,  "ctrlpp::lqr_gain_continuous NX=2",  "ct::optcon::LQR NX=2");
    run_size_sweep<4, 2>(bench,  "ctrlpp::lqr_gain_continuous NX=4",  "ct::optcon::LQR NX=4");
    run_size_sweep<6, 2>(bench,  "ctrlpp::lqr_gain_continuous NX=6",  "ct::optcon::LQR NX=6");
    run_size_sweep<8, 2>(bench,  "ctrlpp::lqr_gain_continuous NX=8",  "ct::optcon::LQR NX=8");
    run_size_sweep<12, 3>(bench, "ctrlpp::lqr_gain_continuous NX=12", "ct::optcon::LQR NX=12");
    run_size_sweep<16, 4>(bench, "ctrlpp::lqr_gain_continuous NX=16", "ct::optcon::LQR NX=16");
    run_size_sweep<20, 5>(bench, "ctrlpp::lqr_gain_continuous NX=20", "ct::optcon::LQR NX=20");
    run_size_sweep<24, 6>(bench, "ctrlpp::lqr_gain_continuous NX=24", "ct::optcon::LQR NX=24");
    run_size_sweep<30, 6>(bench, "ctrlpp::lqr_gain_continuous NX=30", "ct::optcon::LQR NX=30");

    std::ofstream csv("bench_lqr_continuous_vs_ct.csv");
    bench.render(ctrlpp::bench::csv_tpl, csv);
}
