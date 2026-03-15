#include <ctrlpp/pid.h>
#include <ctrlpp/eigen_linalg.h>

#include <cstdio>

int main()
{
    using Policy = ctrlpp::EigenLinalgPolicy;
    using Pid = ctrlpp::Pid<Policy, double, 1, 1, 1>;
    using Vec = Pid::vector_t;

    // Configure a simple PID controller
    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.0);
    cfg.ki = Vec::Constant(1.0);
    cfg.kd = Vec::Constant(0.1);
    cfg.output_min = Vec::Constant(-10.0);
    cfg.output_max = Vec::Constant(10.0);

    Pid pid(cfg);

    // Simulated first-order plant: y(k+1) = a * y(k) + (1-a) * u(k)
    // Time constant ~ 1s, sampled at dt = 0.1s
    constexpr double a = 0.9;
    double y = 0.0;
    constexpr double setpoint = 1.0;
    constexpr double dt = 0.1;
    constexpr int steps = 100;

    std::printf("%-6s  %-8s  %-10s  %-12s  %-10s  %-10s\n",
                "step", "time", "setpoint", "measurement", "control", "error");
    std::printf("------  --------  ----------  ------------  ----------  ----------\n");

    for (int k = 0; k < steps; ++k) {
        Vec sp = Vec::Constant(setpoint);
        Vec meas = Vec::Constant(y);

        auto u = pid.compute(sp, meas, dt);

        double control = u[0];
        double error = setpoint - y;

        std::printf("%-6d  %8.4f  %10.4f  %12.6f  %10.4f  %10.6f\n",
                    k, static_cast<double>(k) * dt, setpoint, y, control, error);

        // Plant update
        y = a * y + (1.0 - a) * control;
    }

    std::printf("\nFinal measurement: %.6f (setpoint: %.1f)\n", y, setpoint);
}
