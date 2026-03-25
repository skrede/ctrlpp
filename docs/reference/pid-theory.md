# PID Control Theory

PID (Proportional-Integral-Derivative) control is the most widely deployed
feedback control algorithm. A PID controller computes a control signal from
the error between a setpoint and a measured process variable, using three
terms that address present error, accumulated past error, and predicted future
error.

## Key Concepts

### Proportional Action

The proportional term produces output proportional to the current error:

```
u_P(t) = Kp * e(t)
```

Increasing Kp reduces steady-state error and speeds up response, but too much
gain causes oscillation and instability.

### Integral Action

The integral term accumulates past error to eliminate steady-state offset:

```
u_I(t) = Ki * integral(e(tau), 0, t)
```

In discrete time this becomes a running sum. Without integral action, most
systems exhibit a constant offset between setpoint and measurement.

### Derivative Action

The derivative term anticipates future error based on its rate of change:

```
u_D(t) = Kd * d/dt e(t)
```

In practice the derivative is never computed on raw error. A first-order
low-pass filter (derivative filter with bandwidth N) suppresses noise
amplification.

### Anti-Windup

When the actuator saturates, the integral term continues accumulating error --
a phenomenon called integrator windup. Anti-windup strategies prevent this:

- **Back-calculation**: feeds the difference between the limited and unlimited
  output back into the integrator with gain Kb
- **Clamping**: stops integration when the output is saturated and the
  integrator would make it worse
- **Conditional integration**: stops integration when the error exceeds a
  threshold

### Tuning Methods

Classical tuning approaches include:

- **Ziegler-Nichols**: based on the ultimate gain and ultimate period,
  provides aggressive initial gains
- **Lambda tuning**: specifies the desired closed-loop time constant
- **Relay autotuning**: uses a relay experiment to estimate the critical point

Modern approaches use optimisation or loop-shaping. In practice, tuning is
iterative and application-specific.

## References

- **Astrom, K. J. and Hagglund, T.** *Advanced PID Control.* ISA, 2006.
  ISBN 978-1-55617-942-6.
  Comprehensive treatment of PID theory, anti-windup strategies, and
  autotuning methods. Primary reference for ctrlpp's PID implementation.

- **Franklin, G. F., Powell, J. D., and Emami-Naeini, A.** *Feedback Control
  of Dynamic Systems.* Pearson, 7th ed., 2015. ISBN 978-0-13-349659-8.
  Standard textbook covering PID within the broader context of classical and
  state-space control.

## Related API Pages

- [PID Controller](../control/pid/README.md) -- policy-based PID with
  compile-time feature selection
- [Anti-Windup Policy](../control/pid/anti-windup.md) -- back-calculation,
  clamping, conditional integration
- [Derivative Filter Policy](../control/pid/derivative-filter.md) --
  first-order derivative low-pass filter
- [PID Composition Guide](../guides/pid/composition.md) -- composing
  multiple policies
- [Your First PID Guide](../guides/intro/your-first-pid.md) -- introductory
  PID tutorial
- [Cascade Control Guide](../guides/pid/cascade.md) -- inner/outer loop
  cascade PID
