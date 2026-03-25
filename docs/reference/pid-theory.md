# PID Control Theory

PID (Proportional-Integral-Derivative) control is the most widely deployed
feedback control algorithm. A PID controller computes a control signal from
the error between a setpoint and a measured process variable, using three
terms that address present error, accumulated past error, and predicted future
error.

## Parallel Form

The continuous-time PID controller in parallel (ISA) form is:

$$
u(t) = K_p \, e(t) + K_i \int_0^t e(\tau) \, d\tau + K_d \frac{de(t)}{dt}
$$

where $e(t) = r(t) - y(t)$ is the error between setpoint $r$ and measurement
$y$, and $K_p$, $K_i$, $K_d$ are the proportional, integral, and derivative
gains respectively.

In the alternative "standard" form the gains are parameterised as:

$$
u(t) = K_p \left[ e(t) + \frac{1}{T_i} \int_0^t e(\tau) \, d\tau + T_d \frac{de(t)}{dt} \right]
$$

where $T_i = K_p / K_i$ is the integral time and $T_d = K_d / K_p$ is the
derivative time.

## Key Concepts

### Proportional Action

The proportional term produces output proportional to the current error:

$$
u_P(t) = K_p \, e(t)
$$

Increasing $K_p$ reduces steady-state error and speeds up response, but too
much gain causes oscillation and instability.

### Integral Action

The integral term accumulates past error to eliminate steady-state offset. In
discrete time with sampling period $T_s$, the backward-Euler approximation is:

$$
u_I[k] = u_I[k-1] + K_i \, T_s \, e[k]
$$

Without integral action, most systems exhibit a constant offset between
setpoint and measurement.

### Derivative Action

The derivative term anticipates future error based on its rate of change. In
practice the derivative is never computed on raw error because differentiation
amplifies high-frequency noise. A first-order low-pass filter with bandwidth
parameter $N$ gives the filtered derivative transfer function:

$$
D(s) = K_d \frac{N s}{s + N}
$$

In discrete time (bilinear/Tustin transform with period $T_s$):

$$
u_D[k] = \frac{2 K_d N}{2 + N T_s} \bigl(e[k] - e[k-1]\bigr) + \frac{2 - N T_s}{2 + N T_s} \, u_D[k-1]
$$

Larger $N$ gives a faster derivative response but less noise attenuation.
Typical values are $N \in [8, 20]$.

### Anti-Windup

When the actuator saturates, the integral term continues accumulating error --
a phenomenon called integrator windup. Anti-windup strategies prevent this:

- **Back-calculation**: feeds the difference between the limited and unlimited
  output back into the integrator with gain $K_b$:

$$
u_I[k] = u_I[k-1] + K_i \, T_s \, e[k] + K_b \bigl(u_{\text{sat}}[k-1] - u[k-1]\bigr)
$$

- **Clamping**: stops integration when the output is saturated and the
  integrator would make it worse, i.e., when
  $\operatorname{sign}(e[k]) = \operatorname{sign}(u[k] - u_{\text{sat}}[k])$

- **Conditional integration**: stops integration when the error exceeds a
  threshold

### Tuning Methods

Classical tuning approaches include:

- **Ziegler-Nichols**: based on the ultimate gain $K_u$ and ultimate period
  $T_u$, provides aggressive initial gains ($K_p = 0.6 K_u$,
  $T_i = 0.5 T_u$, $T_d = 0.125 T_u$)
- **Lambda tuning**: specifies the desired closed-loop time constant
  $\lambda$
- **Relay autotuning**: uses a relay experiment to estimate the critical point

Modern approaches use optimisation or loop-shaping. In practice, tuning is
iterative and application-specific.

## References

- Astrom, K. J. and Hagglund, T. (2006). *Advanced PID Control.* ISA.
  [`astrom2006`]
  Comprehensive treatment of PID theory, anti-windup strategies, and
  autotuning methods. Primary reference for ctrlpp's PID implementation.

- Franklin, G. F., Powell, J. D., and Emami-Naeini, A. (2015). *Feedback
  Control of Dynamic Systems.* Pearson, 7th ed. [`franklin2015`]
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
