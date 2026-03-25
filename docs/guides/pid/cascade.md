# Cascade Control with PID

Cascade control uses an outer (slow) loop whose output becomes the setpoint
for an inner (fast) loop. This structure is common in motion control:
the outer loop tracks position, the inner loop tracks velocity, and only the
inner loop drives the actuator.

## Structure

```
position_sp --> [Outer PID] --> velocity_sp --> [Inner PID] --> torque --> [Plant]
                  ^                               ^                         |
                  |--- position ------------------|--- velocity ------------|
```

The outer loop runs at a lower rate (e.g., 100 Hz) while the inner loop runs
faster (e.g., 1000 Hz). The decimation ratio determines how many inner steps
execute per outer step.

## Complete Program

```cpp
// Usage: ./cascade_pid | gnuplot -p -e "set datafile separator ','; plot '-' skip 1 using 1:3 with lines title 'position', '' using 1:4 with lines title 'velocity'"
#include <ctrlpp/control/pid.h>

#include <iomanip>
#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1, 1, 1>;
    using Vec = Pid::vector_t;

    // Outer loop: position control (slower, larger gains)
    Pid::config_type outer_cfg{};
    outer_cfg.kp = Vec::Constant(5.0);
    outer_cfg.ki = Vec::Constant(1.0);
    outer_cfg.kd = Vec::Constant(0.0);
    outer_cfg.output_min = Vec::Constant(-10.0);
    outer_cfg.output_max = Vec::Constant(10.0);

    // Inner loop: velocity control (faster, smaller gains)
    Pid::config_type inner_cfg{};
    inner_cfg.kp = Vec::Constant(0.1);
    inner_cfg.ki = Vec::Constant(0.5);
    inner_cfg.kd = Vec::Constant(0.0);
    inner_cfg.output_min = Vec::Constant(-1.0);
    inner_cfg.output_max = Vec::Constant(1.0);

    Pid outer(outer_cfg);
    Pid inner(inner_cfg);

    // Plant: simple rotational inertia with damping
    constexpr double J = 0.01;       // inertia
    constexpr double b = 0.1;        // damping
    constexpr double dt_inner = 0.001;
    constexpr double dt_outer = 0.01;
    constexpr int decimation = 10;   // dt_outer / dt_inner
    constexpr double pos_sp = 1.0;
    constexpr double duration = 5.0;

    double position = 0.0;
    double velocity = 0.0;
    double vel_sp = 0.0;
    int counter = 0;

    std::cout << "time,pos_sp,position,velocity,vel_sp,torque\n";

    for (double t = 0.0; t < duration; t += dt_inner)
    {
        // Outer loop executes every decimation steps
        if (counter == 0)
        {
            auto sp = Vec::Constant(pos_sp);
            auto meas = Vec::Constant(position);
            auto tracking = Vec::Constant(velocity);
            auto u_outer = outer.compute(sp, meas, dt_outer, tracking);
            vel_sp = u_outer[0];
        }

        // Inner loop executes every step
        auto sp_inner = Vec::Constant(vel_sp);
        auto meas_inner = Vec::Constant(velocity);
        auto torque_vec = inner.compute(sp_inner, meas_inner, dt_inner);
        double torque = torque_vec[0];

        // Plant dynamics: torque -> acceleration -> velocity -> position
        velocity += (torque - b * velocity) * dt_inner / J;
        position += velocity * dt_inner;

        counter = (counter + 1) % decimation;

        std::cout << std::fixed << std::setprecision(4)
                  << t << "," << pos_sp << "," << position << ","
                  << velocity << "," << vel_sp << "," << torque << "\n";
    }
}
```

## Key Details

- **Decimation** -- the outer loop runs every `decimation` inner steps.
  The inner loop sees a piecewise-constant setpoint between outer updates.

- **Tracking signal** -- `outer.compute(sp, meas, dt, tracking)` accepts an
  optional fourth argument: the actual inner-loop measurement. This enables
  bumpless transfer when the inner loop saturates.

- **Gain separation** -- the outer loop has larger proportional gain to
  command aggressive velocity changes; the inner loop has smaller gains
  appropriate for the faster sampling rate.

## When to Use Cascade

Cascade control is useful when:

- The plant has a fast inner variable (velocity, current) and a slow outer
  variable (position, temperature)
- Disturbances act primarily on the inner variable
- The inner variable is measurable

If only the outer variable is measurable, consider an
[observer-controller composition](../estimation/observer-controller.md)
instead.

## Next Steps

- [PID API Reference](../../control/pid/README.md) -- bumpless transfer and
  tracking mode
- [PID Composition](composition.md) -- adding policies to cascade loops
