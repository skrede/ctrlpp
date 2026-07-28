#include "control_loop_demo.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

int main()
{
    using ctrlpp::control_loop_demo;
    using ctrlpp::kDt;
    using ctrlpp::kSteps;

    const auto demo_opt = control_loop_demo<double>::make();
    if(!demo_opt.has_value())
    {
        std::fprintf(stderr, "lqr_gain refused the plant on host -- %s\n", describe(demo_opt.error()));
        return EXIT_FAILURE;
    }

    const double k0 = (*demo_opt).K(0, 0);
    const double k1 = (*demo_opt).K(0, 1);
    std::printf("gain K = [%.9f, %.9f]\n", k0, k1);

    std::printf("k,t_s,x0,x1,u\n");
    auto demo = *demo_opt;
    for(int k = 0; k <= kSteps; ++k)
    {
        const double x0 = demo.x[0];
        const double x1 = demo.x[1];
        const double u  = demo.step();
        std::printf("%d,%.6f,%.9f,%.9f,%.9f\n", k, k * kDt, x0, x1, u);
    }
    const double final_norm = demo.x.norm();
    std::printf("final_norm = %.9e\n", final_norm);

    // Recommended H753 double tolerance: a safety multiple of the floating-point
    // evaluation-order residual (state update re-associated as B*u + A*x), floored
    // well below the ESP32 float figure so the double leg is held to a tighter bar.
    auto reordered = *demo_opt;
    for(int k = 0; k <= kSteps; ++k)
    {
        const double u  = -(reordered.K * reordered.x)(0);
        reordered.x = reordered.B * u + reordered.A * reordered.x;
    }
    const double residual        = std::fabs(reordered.x.norm() - final_norm);
    const double recommended_tol = std::max(1.0e3 * residual, 1.0e-9);
    std::printf("recommended kH753DoubleTol = %.3e (order residual = %.3e)\n", recommended_tol, residual);

    return EXIT_SUCCESS;
}
