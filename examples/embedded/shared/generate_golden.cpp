#include "golden_verdict.h"
#include "family_analysis.h"
#include "golden_reference.h"
#include "control_loop_demo.h"

#include <cmath>
#include <cstdio>
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
