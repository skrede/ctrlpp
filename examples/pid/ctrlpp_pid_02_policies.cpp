// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./ctrlpp_pid_02_policies' using 1:2 with lines title 'setpoint', '' using 1:3 with lines title 'measurement', '' using 1:4 with lines title 'control'"
// Redirect: ./ctrlpp_pid_02_policies > output.csv

#include "ctrlpp/control/pid.h"

#include <iomanip>
#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1, ctrlpp::anti_windup<ctrlpp::back_calc>, ctrlpp::deriv_filter>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(3.0);
    cfg.ki = Vec::Constant(1.5);
    cfg.kd = Vec::Constant(0.0);
    cfg.output_min = Vec::Constant(-5.0);
    cfg.output_max = Vec::Constant(5.0);
    cfg.template policy<ctrlpp::deriv_filter>().n = {10.0};

    Pid ctrl(cfg);

    double y = 0.0;
    constexpr double setpoint = 1.0;
    constexpr double dt = 0.01;
    constexpr double a = 0.9;
    constexpr double duration = 10.0;

    std::cout << "time,setpoint,measurement,control\n";

    for(double t = 0.0; t < duration; t += dt)
    {
        double disturbance = (t >= 5.0) ? -0.5 : 0.0;

        auto sp = Vec::Constant(setpoint);
        auto meas = Vec::Constant(y + disturbance);
        // A refused cycle produced no command. This example stops, because a
        // refusal here would mean the example itself is wrong. A real caller
        // must instead decide what the actuator does: hold the last command,
        // drive a configured safe value, or fail over.
        auto step = ctrl.compute(sp, meas, dt);
        if(!step.has_value())
        {
            std::cerr << "pid refused the cycle at t=" << t << "\n";
            return 1;
        }
        const auto& u = *step;
        y = a * y + (1.0 - a) * u[0];

        std::cout << std::fixed << std::setprecision(4) << t << "," << setpoint << "," << y << "," << u[0] << "\n";
    }
}
