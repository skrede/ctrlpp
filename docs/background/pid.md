# PID Control

PID (Proportional-Integral-Derivative) control is the most widely deployed
feedback control algorithm in industrial practice. A PID controller computes a
control signal from the error between a desired setpoint and a measured process
variable, combining three terms that address present error, accumulated past
error, and predicted future error [1, Ch. 1, pp. 1--8].

The simplicity and effectiveness of PID control have made it the default choice
for single-loop regulatory control. Despite decades of research into advanced
control strategies, surveys consistently find that over 90% of industrial
control loops use some form of PID [2, Sec. 1.1, pp. 1--3].

## Continuous-Time Formulation

### Parallel (ISA) Form

The continuous-time PID controller in parallel form is [1, eq. (1.1), p. 3]:

$$
u(t) = K_p \, e(t) + K_i \int_0^t e(\tau) \, d\tau + K_d \frac{de(t)}{dt}
$$

where $e(t) = r(t) - y(t)$ is the error between setpoint $r$ and measurement
$y$, and $K_p$, $K_i$, $K_d$ are the proportional, integral, and derivative
gains respectively.

### Standard (Series) Form

An alternative parameterisation uses the proportional gain as a common factor
[1, eq. (1.2), p. 4]:

$$
u(t) = K_p \left[ e(t) + \frac{1}{T_i} \int_0^t e(\tau) \, d\tau + T_d \frac{de(t)}{dt} \right]
$$

where $T_i = K_p / K_i$ is the integral time and $T_d = K_d / K_p$ is the
derivative time. This form makes the physical meaning of each parameter more
transparent: $T_i$ is the time for the integral term to repeat the proportional
action, and $T_d$ is the time by which the derivative term anticipates
the proportional action.

## Proportional Action

The proportional term produces output proportional to the current error
[1, Sec. 1.2.1, p. 5]:

$$
u_P(t) = K_p \, e(t)
$$

Increasing $K_p$ reduces steady-state error and increases the speed of
response, but excessive gain causes oscillation and eventually instability.
For most systems, proportional-only control leaves a residual steady-state
offset because the controller requires a non-zero error to produce a non-zero
output [3, Sec. 9.3, pp. 485--490].

## Integral Action

The integral term accumulates past error, eliminating steady-state offset for
step disturbances and reference changes [1, Sec. 1.2.2, pp. 6--7]. In the
Laplace domain, integral action contributes a pole at the origin, ensuring
infinite DC gain and hence zero steady-state error for constant inputs.

In discrete time with sampling period $T_s$, the backward-Euler approximation
is [1, Sec. 3.2, p. 67]:

$$
u_I[k] = u_I[k-1] + K_i \, T_s \, e[k]
$$

The trapezoidal (Tustin) approximation provides better accuracy [1, Sec. 3.2,
p. 68]:

$$
u_I[k] = u_I[k-1] + \frac{K_i \, T_s}{2} \bigl(e[k] + e[k-1]\bigr)
$$

## Derivative Action

The derivative term anticipates future error based on its rate of change.
In practice, the derivative is never applied to raw error because
differentiation amplifies high-frequency noise. A first-order low-pass filter
with bandwidth parameter $N$ gives the filtered derivative [1, Sec. 1.2.3,
pp. 8--10]:

$$
D(s) = K_d \frac{N s}{s + N}
$$

The parameter $N$ determines the trade-off between derivative response speed
and noise attenuation. Typical values are $N \in [8, 20]$
[1, Sec. 3.3, p. 72].

In discrete time using the bilinear (Tustin) transform with period $T_s$
[1, Sec. 3.3, p. 73]:

$$
u_D[k] = \frac{2 K_d N}{2 + N T_s} \bigl(e[k] - e[k-1]\bigr) + \frac{2 - N T_s}{2 + N T_s} \, u_D[k-1]
$$

## Anti-Windup

When the actuator saturates, the integral term continues accumulating error,
a phenomenon called integrator windup. This leads to large overshoots and
sluggish recovery when the setpoint changes [1, Sec. 3.5, pp. 79--85].

### Back-Calculation

Back-calculation feeds the difference between the limited and unlimited output
back into the integrator with tracking gain $K_b$
[1, Sec. 3.5.1, pp. 80--82]:

$$
u_I[k] = u_I[k-1] + K_i \, T_s \, e[k] + K_b \bigl(u_{\text{sat}}[k-1] - u[k-1]\bigr)
$$

A common choice is $K_b = 1/T_i$, which gives a tracking time constant equal
to the integral time [2, Sec. 3.4, p. 89].

### Clamping (Conditional Integration)

Clamping stops integration when the output is saturated and the integrator
would make the saturation worse [1, Sec. 3.5.2, p. 83]:

$$
\text{if } \operatorname{sign}(e[k]) = \operatorname{sign}(u[k] - u_{\text{sat}}[k]) \implies u_I[k] = u_I[k-1]
$$

This is simpler than back-calculation but can be less smooth during
transitions.

## Setpoint Weighting

Setpoint weighting separates the response to setpoint changes from the
response to disturbances. The proportional and derivative terms use weighted
setpoint signals [1, Sec. 3.6, pp. 86--88]:

$$
e_P(t) = b \cdot r(t) - y(t), \qquad e_D(t) = c \cdot r(t) - y(t)
$$

where $b \in [0, 1]$ and $c \in [0, 1]$ are the proportional and derivative
weights. Setting $b < 1$ reduces overshoot on setpoint steps without
affecting disturbance rejection. Setting $c = 0$ avoids derivative kicks
from step changes in the reference.

## Discrete-Time Implementation

A complete discrete PID implementation combines all the above elements.
The velocity (incremental) form computes the change in output rather than the
absolute output, avoiding the need to store the integral sum explicitly and
providing natural bumpless transfer [1, Sec. 3.4, pp. 75--78]:

$$
\Delta u[k] = u[k] - u[k-1]
$$

This form is preferred in practice because switching between manual and
automatic modes does not cause output bumps.

## Feed-Forward

Pure feedback control reacts only after an error appears. Feed-forward action
uses knowledge of the disturbance or setpoint trajectory to act preemptively
[1, Sec. 5.1, pp. 131--135]:

$$
u(t) = u_{\text{fb}}(t) + u_{\text{ff}}(t)
$$

where $u_{\text{ff}}$ is derived from a model of the process or the known
disturbance path. Feed-forward does not affect stability but can dramatically
improve transient response.

## References

[1] A. Visioli, "Practical PID Control," Springer, 2006.

[2] R. Vilanova and A. Visioli (Eds.), "PID Control in the Third Millennium:
Lessons Learned and New Approaches," Springer, 2012.

[3] N. S. Nise, "Control Systems Engineering," 7th ed., Wiley, 2015.
