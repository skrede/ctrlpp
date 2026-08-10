#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_CT_CT_CARE_ARM_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_CT_CT_CARE_ARM_H

#include "riccati_problem.h"

#include "bench_csv.h"
#include "bench_construct.h"

#include "ctrlpp/control/care.h"

#include <cassert>  // must precede ct_optcon includes; DynamicRiccatiEquation.hpp uses assert() without <cassert>
#include <ct/optcon/lqr/riccati/CARE.hpp>
#include <ct/optcon/lqr/riccati/CARE-impl.hpp>

#include <Eigen/Dense>

#include <cstddef>

namespace ctrlpp::bench
{

constexpr char const* care_deviation_metric = "max abs entrywise deviation of the two arms' Riccati solutions P";
constexpr char const* care_residual_metric = "relative residual of this arm's own continuous Riccati solution";

template <std::size_t NX, std::size_t NU>
class ct_care_arm
{
public:
    explicit ct_care_arm(const riccati_plant<NX, NU>& plant)
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

// The own-criterion rows re-run the identical arm, so their timing columns are a
// second sample of the same work rather than a different workload. Each residual
// is taken on that solver's own solution: routing one solver's P through the
// other's acceptance expression would judge it by a criterion it never met.
template <std::size_t NX, std::size_t NU>
void emit_care_rows(ankerl::nanobench::Bench& bench, const riccati_plant<NX, NU>& plant, ct_care_arm<NX, NU>& ct_arm,
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

    report_accuracy(bench, care_deviation_metric, (P_ctrlpp - P_ct).cwiseAbs().maxCoeff());
    bench.run(label_ctrlpp, solve_ctrlpp).run(label_ct, solve_ct);
    run_own_criterion_pair(bench, care_residual_metric, label_ctrlpp,
                           riccati_relative_residual<NX, NU>(plant, P_ctrlpp), solve_ctrlpp, label_ct,
                           riccati_relative_residual<NX, NU>(plant, P_ct), solve_ct);
}

template <std::size_t NX, std::size_t NU>
void run_care_sweep(ankerl::nanobench::Bench& bench, const char* label_ctrlpp, const char* label_ct)
{
    const riccati_plant<NX, NU> plant = build_damped_chain<NX, NU>();
    const Eigen::Matrix<double, int(NX), int(NX)> P_ctrlpp =
        built_or_exit(ctrlpp::care<double, NX, NU>(plant.A, plant.B, plant.Q, plant.R), label_ctrlpp).P;
    ct_care_arm<NX, NU> ct_arm{plant};
    const Eigen::Matrix<double, int(NX), int(NX)> P_ct = ct_arm.solve();
    emit_care_rows<NX, NU>(bench, plant, ct_arm, label_ctrlpp, label_ct, P_ctrlpp, P_ct);
}

}

#endif
