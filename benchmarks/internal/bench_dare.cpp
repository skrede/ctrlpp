// Discrete Riccati benchmark, swept over state dimension.
//
// The acceptance check is the independent variable. Four quantities are timed
// per size; the first three run on identical operands:
//
//   solve+accept -- the whole `ctrlpp::dare` entry point
//   accept       -- `verify_dare_solution` alone, on the solve's own answer
//   fwd-error    -- `estimate_riccati_forward_error` alone, on pre-formed
//                   closed-loop and residual operands
//   solve+accept scaled
//                -- the same entry point on weights whose scale is s rather
//                   than one, timed at every rung of a geometric ladder in s
//
// The solve without the check is the difference of the first two, which is the
// denominator the published operation counts use: those are stated against a
// solve with no acceptance arithmetic in it. Timing the round's start against
// its end would confound the construction change, the removed factorization and
// the operand preparation into one number that attributes to none of them.
//
// The corpus is the discrete damped chain: a forward-Euler step of the same
// continuous chain the CARE bakeoff sweeps, so the two harnesses describe the
// same family. It is built at two weight scales. Q = I and R = 0.1 I put the
// weight scale at exactly one, so equilibration is inactive and the solve runs
// the acceptance check ONCE. Multiplying both weights by s puts the weight
// scale at exactly s, and the entry point then runs the whole equilibrated
// path: the acceptance check at the equilibrated scale, the rescale of the
// answer, a definiteness factorization of the rescaled answer, the acceptance
// check again at the caller's scale, and a gain agreement between the two
// scales. That is the fourth quantity above.
//
// THE EQUILIBRATED PATH'S COST IS THE DIFFERENCE between the scaled row and the
// scale-one row at the same size. It is not timed separately and it is not
// asserted anywhere, so a reader recomputes it from two published rows rather
// than taking one number on trust.
//
// `--check` runs the corpus confirmation instead of the timing set: it prints,
// for every size and rung, the weight scale the helper computes and whether the
// entry point's equilibrated gate is therefore open. A scaled corpus that
// failed to open that gate would produce a second row equal to the first and an
// overhead of approximately zero, and nothing in a timing would show it.

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include "ctrlpp/detail/riccati_solution.h"
#include "ctrlpp/control/dare.h"

#include <Eigen/Dense>

#include <cstdio>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <tuple>
#include <algorithm>
#include <string_view>

namespace
{

constexpr char const* csv_tpl =
    R"TEMPLATE("title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

/// Forward-Euler step of the continuous damped chain the CARE bakeoff uses:
/// -0.5 on the diagonal and 1.0 on the superdiagonal, one input per group of
/// states. The step keeps every eigenvalue of A at 1 - 0.5 dt, inside the unit
/// disk, so the discrete equation is well posed at every size.
template <std::size_t NX, std::size_t NU>
auto build_discrete_damped_chain()
{
    constexpr int n  = int(NX);
    constexpr int nu = int(NU);
    constexpr double dt = 0.01;

    Eigen::Matrix<double, n, n> A = Eigen::Matrix<double, n, n>::Identity();
    for (std::size_t i = 0; i < NX; ++i)
        A(int(i), int(i)) += dt * -0.5;
    for (std::size_t i = 0; i + 1 < NX; ++i)
        A(int(i), int(i + 1)) = dt * 1.0;

    Eigen::Matrix<double, n, nu> B = Eigen::Matrix<double, n, nu>::Zero();
    const std::size_t group = NX / NU;
    for (std::size_t j = 0; j < NU; ++j)
    {
        const std::size_t last_row = std::min((j + 1) * group, NX) - 1;
        B(int(last_row), int(j)) = dt;
    }

    Eigen::Matrix<double, n, n>   Q = Eigen::Matrix<double, n, n>::Identity();
    Eigen::Matrix<double, nu, nu> R = 0.1 * Eigen::Matrix<double, nu, nu>::Identity();

    return std::tuple{A, B, Q, R};
}

/// The same chain at a weight scale of `s` instead of one.
///
/// `dare_weight_scale` returns max(||Q||_max, ||R||_max), so multiplying both
/// weights by a positive `s` makes it return exactly `s`, and the entry point's
/// gate on that value being different from one opens.
///
/// BOTH WEIGHTS ARE MULTIPLIED AND NOT ONLY Q, deliberately. Scaling Q alone
/// moves the weight ratio with the scale, so the equilibrated problem the solve
/// actually sees would be a different problem at every rung and the difference
/// against the scale-one row would confound the equilibrated path's cost with a
/// change in the work underneath it. With both scaled, Q/s is the identity
/// exactly and R/s reproduces the scale-one corpus to within one rounding of
/// 0.1*s, so the Schur solve underneath is the same work at every rung and the
/// difference is the equilibrated path plus the two rescaling passes.
template <std::size_t NX, std::size_t NU>
auto build_scaled_damped_chain(double s)
{
    auto [A, B, Q, R] = build_discrete_damped_chain<NX, NU>();
    Q *= s;
    R *= s;
    return std::tuple{A, B, Q, R};
}

struct weight_scale_rung
{
    double      scale;
    char const* label;
};

/// A geometric ladder in the weight scale, three decades either side of one,
/// plus one rung immediately beside it.
///
/// The overhead is not scale-independent by construction: dividing the weights
/// by the scale moves the operands' exponents, and both the definiteness
/// factorization's iteration count and the gain agreement's resolution can move
/// with them. A ladder answers whether the cost depends on the scale; one value
/// would answer whether it depends on that value.
///
/// The near-one rung separates two different statements -- "the gate is open"
/// and "the pose is far from equilibrated" -- because everything below the gate
/// runs at 1 + 2^-16 exactly as it runs at 1e+06. That rung is exactly
/// representable, so its weight scale is that value and not a rounding of it.
constexpr weight_scale_rung scale_rungs[] = {
    {1e-06,          "1e-06"  },
    {1e-04,          "1e-04"  },
    {1e-02,          "1e-02"  },
    {1.0 + 0x1p-16,  "1+2^-16"},
    {1e+02,          "1e+02"  },
    {1e+04,          "1e+04"  },
    {1e+06,          "1e+06"  },
};

template <std::size_t NX, std::size_t NU>
void run_size_sweep(ankerl::nanobench::Bench& bench)
{
    constexpr int n  = int(NX);
    constexpr int nu = int(NU);

    auto [A, B, Q, R] = build_discrete_damped_chain<NX, NU>();

    // Solved once, outside every measurement window, so that the acceptance
    // rows are handed the same answer the solve produced rather than one they
    // recompute.
    auto solved = ctrlpp::dare<double, NX, NU>(A, B, Q, R);
    if (!solved)
    {
        std::cerr << "SKIP NX=" << NX << ": dare declined the corpus\n";
        return;
    }
    const Eigen::Matrix<double, n, n>  P = solved->P;
    const Eigen::Matrix<double, nu, n> K = solved->K;

    // The two operands the forward-error estimate consumes, formed exactly as
    // the acceptance path forms them. Building them here rather than inside the
    // timed loop is what makes the fwd-error row the estimate's own cost and
    // not the cost of preparing its inputs.
    const Eigen::Matrix<double, n, n> AtPA        = (A.transpose() * P * A).eval();
    const Eigen::Matrix<double, n, n> AtPBK       = (A.transpose() * P * B * K).eval();
    const Eigen::Matrix<double, n, n> residual    = (AtPA - P - AtPBK + Q).eval();
    const Eigen::Matrix<double, n, n> closed_loop = (A - B * K).eval();

    char buf[64];

    std::snprintf(buf, sizeof(buf), "dare solve+accept NX=%zu", NX);
    bench.run(buf, [&]
    {
        auto r = ctrlpp::dare<double, NX, NU>(A, B, Q, R);
        ankerl::nanobench::doNotOptimizeAway(r);
    });

    std::snprintf(buf, sizeof(buf), "dare accept NX=%zu", NX);
    bench.run(buf, [&]
    {
        Eigen::Matrix<double, nu, n> verified_gain;
        auto v = ctrlpp::detail::verify_dare_solution<double, NX, NU>(A, B, Q, R, P, &verified_gain);
        ankerl::nanobench::doNotOptimizeAway(v);
        ankerl::nanobench::doNotOptimizeAway(verified_gain);
    });

    std::snprintf(buf, sizeof(buf), "dare fwd-error NX=%zu", NX);
    bench.run(buf, [&]
    {
        auto e = ctrlpp::detail::estimate_riccati_forward_error<double, n>(closed_loop, residual, P);
        ankerl::nanobench::doNotOptimizeAway(e);
    });

    // The fourth quantity, one row per rung: the same entry point on weights
    // whose scale is not one, so the call carries the equilibrated path on top
    // of everything the first row already measures. Subtract the two rows at
    // one size to obtain that path's cost.
    for (const weight_scale_rung& rung : scale_rungs)
    {
        auto [A_s, B_s, Q_s, R_s] = build_scaled_damped_chain<NX, NU>(rung.scale);

        std::snprintf(buf, sizeof(buf), "dare solve+accept scaled NX=%zu s=%s", NX, rung.label);
        bench.run(buf, [&]
        {
            auto r = ctrlpp::dare<double, NX, NU>(A_s, B_s, Q_s, R_s);
            ankerl::nanobench::doNotOptimizeAway(r);
        });
    }
}

char const* verification_name(ctrlpp::detail::dare_verification verdict)
{
    switch (verdict)
    {
        case ctrlpp::detail::dare_verification::verified:         return "verified";
        case ctrlpp::detail::dare_verification::refuted:          return "refuted";
        case ctrlpp::detail::dare_verification::unresolved:       return "unresolved";
        case ctrlpp::detail::dare_verification::gain_unavailable: return "gain_unavailable";
    }
    return "unrecognized";
}

/// Confirm, WITHOUT TIMING ANYTHING, that the scaled corpus takes the
/// equilibrated path.
///
/// The gate is one line of `ctrlpp::dare`: it returns as soon as the weight
/// scale is exactly one. So a weight scale different from one, on a call that
/// then succeeds, is the whole of the evidence that the rescale, the
/// definiteness re-check, the caller-scale acceptance check and the gain
/// agreement all ran. Two consequences of that path are printed beside it
/// rather than inferred: the returned solution against the homogeneity relation
/// s * P, which holds only if the rescale executed, and the caller-scale
/// acceptance verdict, which is the check the entry point performs a second
/// time.
template <std::size_t NX, std::size_t NU>
void check_scaled_corpus(std::size_t& rungs_checked, std::size_t& gate_open, std::size_t& failures)
{
    constexpr int n  = int(NX);
    constexpr int nu = int(NU);

    auto [A, B, Q, R] = build_discrete_damped_chain<NX, NU>();

    const double base_scale = ctrlpp::detail::dare_weight_scale<double, NX, NU>(Q, R);
    auto         base       = ctrlpp::dare<double, NX, NU>(A, B, Q, R);
    if (!base)
    {
        std::printf("CHECK NX=%2zu s=base    FAIL: the scale-one corpus was declined\n", NX);
        ++failures;
        return;
    }
    const Eigen::Matrix<double, n, n> P_base = base->P;
    std::printf("CHECK NX=%2zu s=base    weight_scale=%.17g gate=%s\n", NX, base_scale, base_scale != 1.0 ? "open" : "closed");

    for (const weight_scale_rung& rung : scale_rungs)
    {
        ++rungs_checked;

        auto [A_s, B_s, Q_s, R_s] = build_scaled_damped_chain<NX, NU>(rung.scale);

        const double weight_scale = ctrlpp::detail::dare_weight_scale<double, NX, NU>(Q_s, R_s);
        const bool   open         = weight_scale != 1.0;
        if (open)
            ++gate_open;
        else
            ++failures;

        auto scaled = ctrlpp::dare<double, NX, NU>(A_s, B_s, Q_s, R_s);
        if (!scaled)
        {
            std::printf("CHECK NX=%2zu s=%-7s weight_scale=%.17g gate=%-6s FAIL: declined\n", NX, rung.label, weight_scale, open ? "open" : "closed");
            ++failures;
            continue;
        }

        const Eigen::Matrix<double, n, n> P_homogeneous = (rung.scale * P_base).eval();
        const double                      reference     = P_homogeneous.cwiseAbs().maxCoeff();
        const double                      deviation     = (scaled->P - P_homogeneous).cwiseAbs().maxCoeff();
        const double                      relative      = reference > 0.0 ? deviation / reference : deviation;

        Eigen::Matrix<double, nu, n>            caller_gain;
        const ctrlpp::detail::dare_verification caller_verdict =
                ctrlpp::detail::verify_dare_solution<double, NX, NU>(A_s, B_s, Q_s, R_s, scaled->P, &caller_gain);

        std::printf("CHECK NX=%2zu s=%-7s weight_scale=%.17g gate=%-6s rescale_rel_dev=%.3e caller_scale_verdict=%s\n", NX, rung.label, weight_scale,
                    open ? "open" : "closed", relative, verification_name(caller_verdict));
    }
}

/// The size sweep, written once. Both the timing set and the corpus check
/// traverse it, so a confirmation cannot silently cover a different set of
/// sizes than the measurement it licenses.
///
/// NX = 15 IS THE CEILING, AND IT IS A COMPILE-TIME ONE. The forward-error
/// estimate holds an M-by-M operator with M = NX(NX+1)/2 as a fixed-size Eigen
/// object, so Eigen's 128 KiB stack-allocation limit admits M <= 128, that is
/// NX <= 15. At NX = 16 the operator is 136 x 136 = 147,968 bytes and the static
/// assertion fires: the acceptance check cannot be instantiated there at all, so
/// there is no larger size for this sweep to reach.
template <typename Fn>
void for_each_size(Fn&& fn)
{
    fn.template operator()<2, 1>();
    fn.template operator()<4, 2>();
    fn.template operator()<6, 2>();
    fn.template operator()<8, 2>();
    fn.template operator()<12, 3>();
    fn.template operator()<15, 3>();
}

/// AN ARGUMENT THIS PROGRAM DOES NOT RECOGNIZE MUST NOT FALL THROUGH TO THE
/// TIMING SET. The timed rows are comparable to the published figures only when
/// the machine was quiet, and an argument typed by someone expecting a usage
/// message is exactly the invocation that would otherwise produce a timing
/// nobody intended, on a machine nobody had reserved. The measurement runs on no
/// arguments and on nothing else.
int usage(char const* unrecognized)
{
    if (unrecognized != nullptr)
        std::fprintf(stderr, "bench_dare: unrecognized argument '%s'\n", unrecognized);
    std::fprintf(unrecognized != nullptr ? stderr : stdout,
                 "usage: bench_dare           run the timing set (reserved machine only)\n"
                 "       bench_dare --check   confirm the scaled corpus, timing nothing\n"
                 "       bench_dare --help    this message\n");
    return unrecognized != nullptr ? 2 : 0;
}

void check_perf_event_paranoid()
{
    std::ifstream paranoid_file("/proc/sys/kernel/perf_event_paranoid");
    int paranoid_value = 3;
    paranoid_file >> paranoid_value;
    if (paranoid_value > 2)
        std::cerr << "WARN: perf_event_paranoid=" << paranoid_value
                  << "; nanobench hardware counters may be zero. "
                  << "Run:  sudo sysctl kernel.perf_event_paranoid=2\n";
}

}

int main(int argc, char** argv)
{
    bool check_only = false;
    for (int i = 1; i < argc; ++i)
    {
        const std::string_view arg(argv[i]);
        if (arg == "--check")
            check_only = true;
        else if (arg == "--help" || arg == "-h")
            return usage(nullptr);
        else
            return usage(argv[i]);
    }

    // The corpus confirmation, and nothing timed. It belongs BEFORE the
    // measurement rather than inside it: a corpus that does not open the
    // equilibrated gate makes the fourth row a copy of the first, and finding
    // that out during a measurement costs the measurement.
    if (check_only)
    {
        std::size_t rungs_checked = 0;
        std::size_t gate_open     = 0;
        std::size_t failures      = 0;
        for_each_size([&]<std::size_t NX, std::size_t NU>() { check_scaled_corpus<NX, NU>(rungs_checked, gate_open, failures); });
        std::printf("CHECK SUMMARY: rungs=%zu gate_open=%zu failures=%zu\n", rungs_checked, gate_open, failures);
        return failures == 0 ? 0 : 1;
    }

    check_perf_event_paranoid();

    // 51 epochs give an odd-sized sample, so the reported median is an observed
    // measurement rather than an interpolation, and the CSV's error column is
    // the median absolute percent error across those same 51.
    ankerl::nanobench::Bench bench;
    bench.title("DARE acceptance overhead (size sweep)")
        .warmup(10)
        .epochs(51)
        .minEpochTime(std::chrono::milliseconds(1))
        .performanceCounters(true);

    for_each_size([&]<std::size_t NX, std::size_t NU>() { run_size_sweep<NX, NU>(bench); });

    std::ofstream csv("bench_dare.csv");
    bench.render(csv_tpl, csv);
}
