#ifndef HPP_GUARD_BENCHMARKS_INTERNAL_DARE_CORPUS_CHECK_H
#define HPP_GUARD_BENCHMARKS_INTERNAL_DARE_CORPUS_CHECK_H

#include "dare_corpus.h"

#include "ctrlpp/detail/riccati_solution.h"
#include "ctrlpp/control/dare.h"

#include <Eigen/Dense>

#include <cstdio>
#include <cstdint>
#include <cstddef>

namespace ctrlpp::bench
{

inline char const* verification_name(ctrlpp::detail::dare_verification verdict)
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

/// The corpus confirmation, and nothing timed. It belongs BEFORE the
/// measurement rather than inside it: a corpus that does not open the
/// equilibrated gate makes the fourth row a copy of the first, and finding that
/// out during a measurement costs the measurement.
inline int32_t run_corpus_check()
{
    std::size_t rungs_checked = 0;
    std::size_t gate_open     = 0;
    std::size_t failures      = 0;
    for_each_size([&]<std::size_t NX, std::size_t NU>() { check_scaled_corpus<NX, NU>(rungs_checked, gate_open, failures); });
    std::printf("CHECK SUMMARY: rungs=%zu gate_open=%zu failures=%zu\n", rungs_checked, gate_open, failures);
    return failures == 0 ? 0 : 1;
}

}

#endif
