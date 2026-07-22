// ctrlpp NUCLEO-H753ZI on-device control-loop leg (bare superloop, double FPU).
//
// CMSIS startup runs SystemInit + C++ static ctors, then calls main(). This
// drives the shared double control kernel with no RTOS: designs the
// infinite-horizon LQR gain once, runs the closed loop, streams the trajectory
// as CSV over USART3 -> ST-Link VCP, diffs the on-target double result against
// the independent host-double golden, and halts for the operator to read the
// console.

#include "golden_reference.h"
#include "control_loop_demo.h"

#include <Eigen/Dense>

#include <cmath>
#include <cstdio>

namespace ctrlpp
{
void usart3_console_init() noexcept;
}

int main()
{
    ctrlpp::usart3_console_init();

    std::printf("[ctrlpp] NUCLEO-H753ZI bare-superloop control loop (double)\n");

    auto demo = ctrlpp::control_loop_demo<double>::make();
    if(!demo.has_value())
    {
        std::printf("[ctrlpp] lqr_gain FAILED on device -- Riccati did not converge\n");
        for(;;) { }
    }

    const double dev0 = demo->K(0, 0) - ctrlpp::kHostK0;
    const double dev1 = demo->K(0, 1) - ctrlpp::kHostK1;
    std::printf("gain    K = [%.9f, %.9f]\n", demo->K(0, 0), demo->K(0, 1));
    std::printf("host    K = [%.9f, %.9f]\n", ctrlpp::kHostK0, ctrlpp::kHostK1);
    std::printf("double-vs-host gain dev = [%.3e, %.3e]\n", dev0, dev1);

    std::printf("k,t_s,x0,x1,u\n");
    for(int k = 0; k <= ctrlpp::kSteps; ++k)
    {
        const double u = demo->step();
        std::printf("%d,%.4f,%.6f,%.6f,%.6f\n", k, k * ctrlpp::kDt,
                    demo->x[0], demo->x[1], u);
    }

    const double final_norm = demo->x.norm();
    const double err        = std::fabs(final_norm - ctrlpp::kHostFinalNorm);
    std::printf("settling: final |x| = %.3e (host = %.3e, err = %.3e)\n",
                final_norm, ctrlpp::kHostFinalNorm, err);
    std::printf("golden diff %s (tol = %.1e)\n",
                err < ctrlpp::kH753DoubleTol ? "PASS" : "FAIL", ctrlpp::kH753DoubleTol);

    for(;;) { }
    return 0;
}
