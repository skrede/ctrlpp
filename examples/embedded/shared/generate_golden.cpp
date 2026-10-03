#include "golden_verdict.h"
#include "control_loop_demo.h"

#include <cstdio>
#include <cstdlib>

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

}

int main()
{
    auto host   = ctrlpp::control_loop_demo<double>::make();
    auto single = ctrlpp::control_loop_demo<float>::make();
    if(!host.has_value() || !single.has_value())
    {
        std::fprintf(stderr, "lqr_gain refused the plant on the host\n");
        return EXIT_FAILURE;
    }

    std::printf("gain K = [%.9f, %.9f]\n", host->K(0, 0), host->K(0, 1));
    std::printf("k,t_s,x0,x1,u\n");
    run_kernel(*host, true);
    run_kernel(*single, false);
    print_constants(*host);

    ctrlpp::golden_bound bound = ctrlpp::make_golden_bound();
    print_verdict("float", ctrlpp::judge_golden(*single, bound));
    print_verdict("double", ctrlpp::judge_golden(*host, bound));
    return EXIT_SUCCESS;
}
