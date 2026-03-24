// Usage: ./ctrlpp_pid_03_cascade | gnuplot -p -e "plot '-' using 1:4 with lines"
// Redirect: ./ctrlpp_pid_03_cascade > output.csv

#include "ctrlpp/control/pid.h"

#include <iomanip>
#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1, 1, 1>;
    using Vec = Pid::vector_t;

    Pid::config_type outer_cfg{};
    outer_cfg.kp = Vec::Constant(5.0);
    outer_cfg.ki = Vec::Constant(1.0);
    outer_cfg.kd = Vec::Constant(0.0);
    outer_cfg.output_min = Vec::Constant(-10.0);
    outer_cfg.output_max = Vec::Constant(10.0);

    Pid::config_type inner_cfg{};
    inner_cfg.kp = Vec::Constant(0.1);
    inner_cfg.ki = Vec::Constant(0.5);
    inner_cfg.kd = Vec::Constant(0.0);
    inner_cfg.output_min = Vec::Constant(-1.0);
    inner_cfg.output_max = Vec::Constant(1.0);

    Pid outer(outer_cfg);
    Pid inner(inner_cfg);

    constexpr double J = 0.01;
    constexpr double b = 0.1;
    constexpr double dt_inner = 0.001;
    constexpr double dt_outer = 0.01;
    constexpr int decimation = 10;
    constexpr double pos_sp = 1.0;
    constexpr double duration = 5.0;

    double position = 0.0;
    double velocity = 0.0;
    double vel_sp = 0.0;
    int counter = 0;

    std::cout << "time,pos_sp,position,velocity,vel_sp,torque\n";

    for(double t = 0.0; t < duration; t += dt_inner)
    {
        if(counter == 0)
        {
            auto sp = Vec::Constant(pos_sp);
            auto meas = Vec::Constant(position);
            auto tracking = Vec::Constant(velocity);
            auto u_outer = outer.compute(sp, meas, dt_outer, tracking);
            vel_sp = u_outer[0];
        }

        auto sp_inner = Vec::Constant(vel_sp);
        auto meas_inner = Vec::Constant(velocity);
        auto torque_vec = inner.compute(sp_inner, meas_inner, dt_inner);
        double torque = torque_vec[0];

        velocity += (torque - b * velocity) * dt_inner / J;
        position += velocity * dt_inner;

        counter = (counter + 1) % decimation;

        std::cout << std::fixed << std::setprecision(4) << t << "," << pos_sp << "," << position << "," << velocity << "," << vel_sp << "," << torque << "\n";
    }
}
