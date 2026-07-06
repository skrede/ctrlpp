// pid_linear_step.cpp -- Linear PI controller output matching pid_linear_step.m
// Open-loop: feeds known error signal to PID, compares control output.
// Usage: ./pid_linear_step > pid_linear_step_cpp.csv

#include "ctrlpp/control/pid.h"

#include <cmath>
#include <cstdio>
#include <numbers>

int main()
{
    using Pid = ctrlpp::pid<double, 1>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.0);
    cfg.ki = Vec::Constant(1.0);
    cfg.kd = Vec::Constant(0.0);

    Pid ctrl(cfg);

    constexpr double dt = 0.01;
    constexpr int n_steps = 500;

    std::printf("time,error,control\n");

    for(int k = 0; k < n_steps; ++k)
    {
        double t = k * dt;
        double e = std::exp(-0.5 * t) * std::sin(2.0 * std::numbers::pi * 0.5 * t);

        // PID compute: setpoint = e, measurement = 0  =>  error = e
        auto sp = Vec::Constant(e);
        auto meas = Vec::Constant(0.0);
        auto u = ctrl.compute(sp, meas, dt);

        std::printf("%.15e,%.15e,%.15e\n", t, e, u[0]);
    }
}
