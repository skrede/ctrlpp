// Usage: ./ctrlpp_basic_pid_console | gnuplot -p -e "plot '-' using 1:3 with lines"
// Redirect: ./ctrlpp_basic_pid_console > output.csv

#include "ctrlpp/pid.h"

#include <iomanip>
#include <iostream>

int main()
{
    using Pid = ctrlpp::Pid<double, 1, 1, 1>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.0);
    cfg.ki = Vec::Constant(1.0);
    cfg.kd = Vec::Constant(0.1);
    cfg.output_min = Vec::Constant(-10.0);
    cfg.output_max = Vec::Constant(10.0);

    Pid pid(cfg);

    double y = 0.0;
    constexpr double setpoint = 1.0;
    constexpr double dt = 0.01;
    constexpr double a = 0.9;
    constexpr double duration = 10.0;

    std::cout << "time,setpoint,measurement,control" << std::endl;

    for (double t = 0.0; t < duration; t += dt) {
        Vec sp = Vec::Constant(setpoint);
        Vec meas = Vec::Constant(y);
        auto u = pid.compute(sp, meas, dt);
        y = a * y + (1.0 - a) * u[0];

        std::cout << std::fixed << std::setprecision(4)
                  << t << "," << setpoint << "," << y << "," << u[0] << std::endl;
    }
}
