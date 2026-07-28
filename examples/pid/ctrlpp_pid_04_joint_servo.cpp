// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./ctrlpp_pid_04_joint_servo' using 1:2 with lines title 'pos setpoint', '' using 1:3 with lines title 'position', '' using 1:6 with lines title 'torque'"
// Redirect: ./ctrlpp_pid_04_joint_servo > output.csv

#include "ctrlpp/control/pid.h"

#include <iomanip>
#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1>;
    using Vec = Pid::vector_t;

    Pid::config_type outer_cfg{};
    outer_cfg.kp = Vec::Constant(6.0);
    outer_cfg.ki = Vec::Constant(0.5);
    outer_cfg.kd = Vec::Constant(0.0);
    outer_cfg.output_min = Vec::Constant(-5.0);
    outer_cfg.output_max = Vec::Constant(5.0);

    Pid::config_type inner_cfg{};
    inner_cfg.kp = Vec::Constant(0.01);
    inner_cfg.ki = Vec::Constant(0.1);
    inner_cfg.kd = Vec::Constant(0.0);
    inner_cfg.output_min = Vec::Constant(-2.0);
    inner_cfg.output_max = Vec::Constant(2.0);

    Pid outer(outer_cfg);
    Pid inner(inner_cfg);

    constexpr double J = 0.1;
    constexpr double b = 0.5;
    constexpr double gear_ratio = 50.0;
    constexpr double coulomb_friction = 0.02;
    constexpr double dt_inner = 0.001;
    constexpr double dt_outer = 0.01;
    constexpr int decimation = 10;
    constexpr double target_pos = 1.5707963267948966;
    constexpr double ramp_time = 2.0;
    constexpr double duration = 5.0;

    double position = 0.0;
    double velocity = 0.0;
    double vel_sp = 0.0;
    int counter = 0;

    std::cout << "time,pos_sp,position,velocity,vel_sp,torque\n";

    for(double t = 0.0; t < duration; t += dt_inner)
    {
        double pos_sp = (t < ramp_time) ? target_pos * t / ramp_time : target_pos;

        if(counter == 0)
        {
            auto sp = Vec::Constant(pos_sp);
            auto meas = Vec::Constant(position);
            auto tracking = Vec::Constant(velocity);
            auto u_outer = outer.compute(sp, meas, dt_outer, tracking);
            if(!u_outer.has_value())
            {
                std::cerr << "outer loop refused the cycle at t=" << t << "\n";
                return 1;
            }
            vel_sp = (*u_outer)[0];
        }

        auto sp_inner = Vec::Constant(vel_sp);
        auto meas_inner = Vec::Constant(velocity);
        // A refused cycle produced no command. This example stops, because a
        // refusal here would mean the example itself is wrong. A real caller
        // must instead decide what the actuator does: hold the last command,
        // drive a configured safe value, or fail over. In a cascade the outer
        // loop's refusal also leaves the inner loop without a fresh setpoint,
        // so the decision has to cover both.
        auto torque_vec = inner.compute(sp_inner, meas_inner, dt_inner);
        if(!torque_vec.has_value())
        {
            std::cerr << "inner loop refused the cycle at t=" << t << "\n";
            return 1;
        }
        double torque = (*torque_vec)[0];

        double friction = coulomb_friction * ((velocity > 0.0) ? 1.0 : (velocity < 0.0) ? -1.0 : 0.0);
        velocity += (gear_ratio * torque - b * velocity - friction) * dt_inner / J;
        position += velocity * dt_inner;

        counter = (counter + 1) % decimation;

        std::cout << std::fixed << std::setprecision(4) << t << "," << pos_sp << "," << position << "," << velocity << "," << vel_sp << "," << torque << "\n";
    }
}
