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
#include "ctrlpp/pid.h"

#include <iomanip>
#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1, 1, 1>;
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
        auto u = ctrl.compute(sp, meas, dt);

        // Plant dynamics
        y = a * y + (1.0 - a) * u[0];

        std::cout << std::fixed << std::setprecision(4)
                  << t << "," << setpoint << "," << y << "," << u[0] << "\n";
    }
}
```

## What Is Happening

1. **Template parameters** `pid<double, 1, 1, 1>` -- scalar type `double`,
   one state, one input, one output. This is the SISO specialisation.

2. **Configuration** -- `kp`, `ki`, `kd` are `Eigen::Vector` types (here
   1-dimensional). Output limits prevent actuator saturation.

3. **Control loop** -- `ctrl.compute(setpoint, measurement, dt)` returns the
   control signal. The plant model advances one step, and the loop repeats.

4. **CSV output** -- pipe to gnuplot or load in a spreadsheet to visualise the
   step response.

## Adding a Policy

The bare `pid` has no anti-windup or derivative filtering. Adding a policy is
a compile-time template parameter -- zero runtime cost when not used:

```cpp
using Pid = ctrlpp::pid<double, 1, 1, 1,
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
[PID Composition](../pid/composition.md) for a deep dive on composing
multiple policies together.

## Next Steps

- [PID API Reference](../../control/pid/README.md) -- full method signatures
  and all config fields
- [PID Composition Guide](../pid/composition.md) -- composing anti-windup,
  derivative filtering, rate limiting and more
- [PID Theory](../../reference/pid-theory.md) -- the math behind PID control
