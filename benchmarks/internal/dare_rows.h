#ifndef HPP_GUARD_BENCHMARKS_INTERNAL_DARE_ROWS_H
#define HPP_GUARD_BENCHMARKS_INTERNAL_DARE_ROWS_H

#include "dare_corpus.h"

#include "bench_csv.h"

#include "comparison/ct/riccati_problem.h"

#include "ctrlpp/detail/riccati_solution.h"
#include "ctrlpp/control/dare.h"

#include <Eigen/Dense>

#include <cstdio>
#include <cstddef>
#include <iostream>

namespace ctrlpp::bench
{

/// The two operands the forward-error estimate consumes, formed exactly as the
/// acceptance path forms them. Building them outside the timed loop is what
/// makes the fwd-error row the estimate's own cost and not the cost of
/// preparing its inputs.
template <std::size_t NX>
struct forward_error_operands
{
    Eigen::Matrix<double, int(NX), int(NX)> residual;
    Eigen::Matrix<double, int(NX), int(NX)> closed_loop;
};

template <std::size_t NX, std::size_t NU>
forward_error_operands<NX> build_forward_error_operands(const ctrlpp::bench::riccati_plant<NX, NU>& plant,
                                                        const Eigen::Matrix<double, int(NX), int(NX)>& P,
                                                        const Eigen::Matrix<double, int(NU), int(NX)>& K)
{
    constexpr int n = int(NX);
    const Eigen::Matrix<double, n, n> AtPA  = (plant.A.transpose() * P * plant.A).eval();
    const Eigen::Matrix<double, n, n> AtPBK = (plant.A.transpose() * P * plant.B * K).eval();
    return {(AtPA - P - AtPBK + plant.Q).eval(), (plant.A - plant.B * K).eval()};
}

/// The fourth quantity, one row per rung: the same entry point on weights whose
/// scale is not one, so the call carries the equilibrated path on top of
/// everything the first row already measures. Subtract the scaled row from the
/// scale-one row at the same size to obtain that path's cost.
template <std::size_t NX, std::size_t NU>
void run_scaled_rows(ankerl::nanobench::Bench& bench)
{
    char buf[64];
    for (const weight_scale_rung& rung : scale_rungs)
    {
        auto [A_s, B_s, Q_s, R_s] = build_scaled_damped_chain<NX, NU>(rung.scale);
        auto solved = ctrlpp::dare<double, NX, NU>(A_s, B_s, Q_s, R_s);
        if (!solved)
        {
            std::cerr << "SKIP NX=" << NX << " s=" << rung.label << ": dare declined the scaled corpus\n";
            continue;
        }
        auto solve = [&] {
            auto r = ctrlpp::dare<double, NX, NU>(A_s, B_s, Q_s, R_s);
            ankerl::nanobench::doNotOptimizeAway(r);
        };
        std::snprintf(buf, sizeof(buf), "dare solve+accept scaled NX=%zu s=%s", NX, rung.label);
        ctrlpp::bench::run_single_implementation_row(bench, buf, solve);
        ctrlpp::bench::run_own_criterion_row(
            bench, dare_residual_metric, buf,
            ctrlpp::bench::dare_relative_residual<NX, NU>(
                ctrlpp::bench::riccati_plant<NX, NU>{A_s, B_s, Q_s, R_s}, solved->P),
            solve);
    }
}

/// The acceptance check and the forward-error estimate answer a verdict and a
/// bound rather than a Riccati equation, so neither has a residual of its own.
template <std::size_t NX, std::size_t NU>
void run_acceptance_rows(ankerl::nanobench::Bench& bench, const ctrlpp::bench::riccati_plant<NX, NU>& plant,
                         const Eigen::Matrix<double, int(NX), int(NX)>& P,
                         const forward_error_operands<NX>& operands)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "dare accept NX=%zu", NX);
    ctrlpp::bench::run_single_implementation_row(bench, buf, [&]
    {
        Eigen::Matrix<double, int(NU), int(NX)> verified_gain;
        auto v = ctrlpp::detail::verify_dare_solution<double, NX, NU>(plant.A, plant.B, plant.Q, plant.R, P,
                                                                      &verified_gain);
        ankerl::nanobench::doNotOptimizeAway(v);
        ankerl::nanobench::doNotOptimizeAway(verified_gain);
    });

    std::snprintf(buf, sizeof(buf), "dare fwd-error NX=%zu", NX);
    ctrlpp::bench::run_single_implementation_row(bench, buf, [&]
    {
        auto e = ctrlpp::detail::estimate_riccati_forward_error<double, int(NX)>(operands.closed_loop,
                                                                                operands.residual, P);
        ankerl::nanobench::doNotOptimizeAway(e);
    });
}

template <std::size_t NX, std::size_t NU>
void run_size_sweep(ankerl::nanobench::Bench& bench)
{
    auto [A, B, Q, R] = build_discrete_damped_chain<NX, NU>();
    const ctrlpp::bench::riccati_plant<NX, NU> plant{A, B, Q, R};

    // Solved once, outside every measurement window, so that the acceptance
    // rows are handed the same answer the solve produced rather than one they
    // recompute.
    auto solved = ctrlpp::dare<double, NX, NU>(A, B, Q, R);
    if (!solved)
    {
        std::cerr << "SKIP NX=" << NX << ": dare declined the corpus\n";
        return;
    }
    const Eigen::Matrix<double, int(NX), int(NX)> P = solved->P;

    char buf[64];
    auto solve = [&] {
        auto r = ctrlpp::dare<double, NX, NU>(A, B, Q, R);
        ankerl::nanobench::doNotOptimizeAway(r);
    };
    std::snprintf(buf, sizeof(buf), "dare solve+accept NX=%zu", NX);
    ctrlpp::bench::run_single_implementation_row(bench, buf, solve);
    ctrlpp::bench::run_own_criterion_row(bench, dare_residual_metric, buf,
                                         ctrlpp::bench::dare_relative_residual<NX, NU>(plant, P), solve);

    run_acceptance_rows<NX, NU>(bench, plant, P,
                                build_forward_error_operands<NX, NU>(plant, P, solved->K));
    run_scaled_rows<NX, NU>(bench);
}

}

#endif
