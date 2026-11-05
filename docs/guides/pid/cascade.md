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
    using Pid = ctrlpp::pid<double, 1>;
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

    auto outer_created = Pid::create(outer_cfg);
    auto inner_created = Pid::create(inner_cfg);
    if (!outer_created.has_value() || !inner_created.has_value()) return 1;
    auto& outer = *outer_created;
    auto& inner = *inner_created;

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
            // While the inner loop is saturated, the outer loop tracks the
            // velocity actually reached instead of winding up.
            auto tracking = Vec::Constant(velocity);
            auto u_outer = inner.saturated()
                    ? outer.compute(sp, meas, dt_outer, tracking)
                    : outer.compute(sp, meas, dt_outer);
            if (!u_outer.has_value()) return 1;
            vel_sp = (*u_outer)[0];
        }

        // Inner loop executes every step
        auto sp_inner = Vec::Constant(vel_sp);
        auto meas_inner = Vec::Constant(velocity);
        // A refused cycle produced no command; the caller decides what the
        // actuator does. In a cascade the outer loop's refusal also leaves
        // the inner loop without a fresh setpoint, so the decision covers both.
        auto torque_vec = inner.compute(sp_inner, meas_inner, dt_inner);
        if (!torque_vec.has_value()) return 1;
        double torque = (*torque_vec)[0];

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

- **Decimation**<br/>the outer loop runs every `decimation` inner steps.
  The inner loop sees a piecewise-constant setpoint between outer updates.

- **Tracking signal**<br/>the overload `outer.compute(sp, meas, dt, tracking)`
  re-seeds the outer integrator after the cycle to the tracking signal
  (clamped to the outer limits) less the cycle's proportional and derivative
  terms, the value with which that cycle's command would have equaled the
  signal. It does this on every call and does not
  look at the inner loop itself, so the program calls it only while
  `inner.saturated()` reports that the inner loop could not follow its
  setpoint, with the velocity actually reached as the signal. The outer
  integrator then stays consistent with what the inner loop achieves instead
  of winding up. Calling it on every cycle would replace the outer integrator
  each time and leave the outer loop without integral action of its own. The
  tracking signal is checked before the cycle runs, because it is
  back-assigned into the integrator once the cycle succeeds; a non-finite one
  is rejected as `pid_step_error::non_finite_tracking_signal` and the
  integrator is left untouched.

- **Gain separation**<br/>the outer gains map position error to a velocity
  command, the inner gains velocity error to torque, so their sizes are not
  comparable. The integral is formed as `ki * e * dt`, so the gains are
  continuous-time values and do not change with the sampling rate.

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

- [PID API Reference](../../api/control/pid/README.md)<br/> bumpless transfer and
  tracking mode
- [PID Composition](composition.md)<br/> adding policies to cascade loops
