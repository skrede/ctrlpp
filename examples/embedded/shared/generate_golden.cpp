#include "golden_verdict.h"
#include "family_analysis.h"
#include "golden_reference.h"
#include "control_loop_demo.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <cstdlib>
#include <cstddef>

namespace {

template<class Scalar>
void run_kernel(ctrlpp::control_loop_demo<Scalar> &demo, bool print_rows)
{
    for(int k = 0; k <= ctrlpp::kSteps; ++k)
    {
        const double x0 = static_cast<double>(demo.x[0]);
        const double x1 = static_cast<double>(demo.x[1]);
        const double u  = static_cast<double>(demo.step());
        if(print_rows)
            std::printf("%d,%.6f,%.9f,%.9f,%.9f\n", k, k * ctrlpp::kDt, x0, x1, u);
    }
}

void print_constants(const ctrlpp::control_loop_demo<double> &host)
{
    std::printf("kHostK0        = %.17g\n", host.K(0, 0));
    std::printf("kHostK1        = %.17g\n", host.K(0, 1));
    std::printf("kHostFinalNorm = %.17g\n", static_cast<double>(host.x.norm()));
}

// The observation is printed beneath each bound and never feeds it; the gain
// departure the trajectory bound consumes is gated by its own counted bound.
void print_verdict(const char *leg, const ctrlpp::golden_verdict &verdict)
{
    std::printf("%s gain departure = [%.3e, %.3e], bound = %.3e\n", leg, verdict.gain_departure[0], verdict.gain_departure[1], verdict.gain_tolerance);
    std::printf("%s final-norm departure = %.3e, bound = %.3e (%.3e of the reference norm)\n", leg, verdict.norm_departure, verdict.norm_tolerance,
                verdict.norm_tolerance / ctrlpp::kHostFinalNorm);
    std::printf("%s verdict = %s\n", leg, verdict.pass ? "PASS" : "FAIL");
}

// The bounds validation/cases/lqr_closed_loop_settling/tolerance.cfg carries,
// for ctrlpp's arm only; Octave's arm is not counted. Each gain entry is within
// the counted gain bound of the exact design, so the gain row departs by at most
// nx times it, and the run is bounded against the exact loop, whose norms are
// evaluated at the host's gain. Every column also carries both arms' printing at
// kPrintedDecimals digits after the point and the comparator's reading back.
void print_cross_validation_tolerances()
{
    constexpr int kPrintedDecimals = 15;
    const double eps               = std::numeric_limits<double>::epsilon();
    const double gain_norm         = std::abs(ctrlpp::kHostK0) + std::abs(ctrlpp::kHostK1);
    const double gain              = ctrlpp::gain_departure_bound<double>(ctrlpp::kGoldenNx, ctrlpp::kGoldenNu, gain_norm);
    ctrlpp::golden_bound bound     = ctrlpp::make_golden_bound();
    const double state             = bound.final_state_departure(eps, ctrlpp::kPlantRoundings, static_cast<double>(ctrlpp::kGoldenNx) * gain);
    const double printing          = 2.0 * (0.5 * std::pow(10.0, -kPrintedDecimals) + eps / 2.0);
    const double norm_roundings    = static_cast<double>(ctrlpp::kGoldenNx + 1) * eps;
    std::printf("cross-validation: gain_norm = %.17g, gain bound = %.17g, state bound = %.17g\n", gain_norm, gain, state);
    std::printf("atol_K_00=%.17g\nrtol_K_00=%.17g\natol_K_01=%.17g\nrtol_K_01=%.17g\n", gain, printing, gain, printing);
    std::printf("atol_x0_final=%.17g\nrtol_x0_final=%.17g\natol_x1_final=%.17g\nrtol_x1_final=%.17g\n", state, printing, state, printing);
    std::printf("atol_final_norm=%.17g\nrtol_final_norm=%.17g\n", std::sqrt(static_cast<double>(ctrlpp::kGoldenNx)) * state, printing + norm_roundings);
}

bool report_control()
{
    auto host   = ctrlpp::control_loop_demo<double>::make();
    auto single = ctrlpp::control_loop_demo<float>::make();
    if(!host.has_value() || !single.has_value())
    {
        std::fprintf(stderr, "lqr_gain refused the plant on the host\n");
        return false;
    }

    std::printf("gain K = [%.9f, %.9f]\n", host->K(0, 0), host->K(0, 1));
    std::printf("k,t_s,x0,x1,u\n");
    run_kernel(*host, true);
    run_kernel(*single, false);
    print_constants(*host);

    ctrlpp::golden_bound bound = ctrlpp::make_golden_bound();
    print_verdict("float", ctrlpp::judge_golden(*single, bound));
    print_verdict("double", ctrlpp::judge_golden(*host, bound));
    print_cross_validation_tolerances();
    return true;
}

struct family_row
{
    const char *name;
    const char *constant;
    const char *scale_constant;
    ctrlpp::family_analysis analysis;
    double single;
    double roundings;
    double header_reference;
};

// The generator fails on a family whose count is not a bound for this run,
// whose header reference lies further from this run than two double runs may
// part, or whose float run falls outside its own bound.
bool report_family(const family_row &row)
{
    const ctrlpp::family_analysis &a = row.analysis;
    const double departure           = std::abs(row.single - a.reference);
    const double bound               = ctrlpp::paired_departure_bound<float>(row.roundings, a.scale);
    const bool premise               = a.premise <= a.premise_limit;
    const double agreement           = ctrlpp::paired_departure_bound<double>(row.roundings, a.scale);
    const bool current               = ctrlpp::within_bound(std::abs(a.reference - row.header_reference), agreement);
    std::printf("%s = %.17g\n%s = %.17g\n", row.constant, a.reference, row.scale_constant, a.scale);
    std::printf("%s factors: roundings = %.0f, scale = %.6e, %s = %.17g\n", row.name, row.roundings, a.scale, a.measured_name, a.measured);
    std::printf("%s premise: %s = %.6g, at most %.6g: %s\n", row.name, a.premise_name, a.premise, a.premise_limit, premise ? "holds" : "FAILS");
    std::printf("%s float departure = %.3e, bound = %.3e (%.3e of the reference); double bound = %.3e\n", row.name, departure, bound, bound / std::abs(a.reference), agreement);
    std::printf("%s header reference %s, verdict = %s\n", row.name, current ? "agrees within the double bound" : "is STALE", ctrlpp::within_bound(departure, bound) ? "PASS" : "FAIL");
    return premise && current && ctrlpp::within_bound(departure, bound);
}

template<class Demo>
double run_to_end(Demo demo)
{
    double last = 0.0;
    for(std::size_t k = 0; k < ctrlpp::kGoldenSteps; ++k)
        last = static_cast<double>(demo.step());
    return last;
}

template<template<class> class Demo, class Analyze>
bool report(family_row row, Analyze analyze)
{
    auto host   = Demo<double>::make();
    auto single = Demo<float>::make();
    if(!host.has_value() || !single.has_value())
    {
        std::fprintf(stderr, "%s construction refused on the host\n", row.name);
        return false;
    }
    row.analysis = analyze(*host);
    row.single   = run_to_end(*single);
    return report_family(row);
}

}

int main()
{
    using namespace ctrlpp;
    const bool control = report_control();
    const bool estimation =
            report<estimation_demo>({"estimation", "kHostKalmanVelocityVariance", "kHostKalmanScale", {}, 0.0, estimation_roundings(), kHostKalmanVelocityVariance}, analyze_estimation);
    const bool dsp = report<dsp_demo>({"dsp", "kHostBiquadOutput", "kHostBiquadScale", {}, 0.0, dsp_roundings(), kHostBiquadOutput}, analyze_dsp);
    const bool trajectory =
            report<trajectory_demo>({"trajectory", "kHostTrajectoryPosition", "kHostTrajectoryScale", {}, 0.0, trajectory_roundings(), kHostTrajectoryPosition}, analyze_trajectory);
    return control && estimation && dsp && trajectory ? EXIT_SUCCESS : EXIT_FAILURE;
}
