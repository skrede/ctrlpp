# Your First PID Controller

This tutorial walks through building a complete PID control loop from scratch.
By the end, you will have a working controller that regulates a simple
first-order plant (a DC motor model) to a desired setpoint.

**Prerequisites:** ctrlpp installed per [Getting Started](../../getting-started.md).

## The Plant

We model a DC motor as a discrete first-order lag. At each time step the
output blends the previous output with the control signal:

```
y[k+1] = a * y[k] + (1 - a) * u[k]
```

where `a = 0.9` sets the time constant.

## Complete Program

```cpp
// Usage: ./your_first_pid | gnuplot -p -e "set datafile separator ','; plot '-' skip 1 using 1:3 with lines title 'measurement', '' using 1:2 with lines title 'setpoint'"
#include <ctrlpp/control/pid.h>

#include <iomanip>
#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1>;
    using Vec = Pid::vector_t;

    // Configure gains and output limits
    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.0);
    cfg.ki = Vec::Constant(1.0);
    cfg.kd = Vec::Constant(0.0);
    cfg.output_min = Vec::Constant(-10.0);
    cfg.output_max = Vec::Constant(10.0);

    Pid ctrl(cfg);

    // Simulated first-order plant
    double y = 0.0;
    constexpr double setpoint = 1.0;
    constexpr double dt = 0.01;
    constexpr double a = 0.9;
    constexpr double duration = 10.0;

    std::cout << "time,setpoint,measurement,control\n";

    for (double t = 0.0; t < duration; t += dt)
    {
        auto sp = Vec::Constant(setpoint);
        auto meas = Vec::Constant(y);
        // A refused cycle produced no command; the caller decides what the
        // actuator does. This sample stops.
        auto step = ctrl.compute(sp, meas, dt);
        if (!step.has_value()) return 1;
        const auto& u = *step;

        // Plant dynamics
        y = a * y + (1.0 - a) * u[0];

        std::cout << std::fixed << std::setprecision(4)
                  << t << "," << setpoint << "," << y << "," << u[0] << "\n";
    }
}
```

## What Is Happening

1. **Template parameters** `pid<double, 1>`<br/> Scalar type `double`,
   one state, one input, one output. This is the SISO specialization.

2. **Configuration**<br/>`kp`, `ki`, `kd` are `Eigen::Vector` types (here
   1-dimensional). Output limits prevent actuator saturation.

3. **Control loop**<br/>`ctrl.compute(setpoint, measurement, dt)` returns a
   result carrying the control signal, because a cycle can fail: a non-finite
   setpoint or measurement, or a step that is not a positive finite duration, is
   rejected before any carried state is touched. A refusal means there is no
   command for this cycle, so the caller decides what the actuator does -- hold
   the last successful command, drive a configured safe value, or fail over. On
   success the plant model advances one step and the loop repeats.

4. **CSV output**<br/>Pipe to gnuplot or load in a spreadsheet to visualise the
   step response.

## Adding a Policy

The bare `pid` has no anti-windup or derivative filtering. Adding a policy is
a compile-time template parameter &mdash; zero runtime cost when not used:

```cpp
#include <ctrlpp/control/pid.h>

using Pid = ctrlpp::pid<double, 1,
                         ctrlpp::anti_windup<ctrlpp::back_calc>,
                         ctrlpp::deriv_filter>;

Pid::config_type cfg{};
cfg.kp = Vec::Constant(3.0);
cfg.ki = Vec::Constant(1.5);
cfg.kd = Vec::Constant(0.5);
cfg.output_min = Vec::Constant(-5.0);
cfg.output_max = Vec::Constant(5.0);

// Configure the derivative filter bandwidth
cfg.template policy<ctrlpp::deriv_filter>().n = {10.0};
```

The config struct automatically gains fields for each policy. See
[PID Composition](../pid/composition.md) for a details on composing
multiple policies together.

## Next Steps

- [PID API Reference](../../api/control/pid/README.md)<br/> full method signatures
  and all config fields
- [PID Composition Guide](../pid/composition.md)<br/> composing anti-windup,
  derivative filtering, rate limiting and more
- [PID Theory](../../background/pid.md)<br/> the math behind PID control
