# PID Control

PID (Proportional-Integral-Derivative) control is the most widely deployed
feedback control algorithm in industrial practice. A PID controller computes a
control signal from the error between a desired setpoint and a measured process
variable, combining three terms that address present error, accumulated past
error, and predicted future error [1, Ch. 1, pp. 1--8].

The simplicity and effectiveness of PID control have made it the default choice
for single-loop regulatory control. Despite decades of research into advanced
control strategies, surveys consistently find that over 90% of industrial
control loops use some form of PID [2, Sec. 1.1, p. 3].

## Continuous-Time Formulation

### Parallel Form

The continuous-time PID controller in parallel form is [1, eq. (1.18), p. 8]:

$$
u(t) = K_p \, e(t) + K_i \int_0^t e(\tau) \, d\tau + K_d \frac{de(t)}{dt}
$$

where $e(t) = r(t) - y(t)$ is the error between setpoint $r$ and measurement
$y$, and $K_p$, $K_i$, $K_d$ are the proportional, integral, and derivative
gains respectively. This is the form `pid` computes in, with the gains in
`kp`, `ki` and `kd`.

### Ideal (ISA) Form

An alternative parametrization uses the proportional gain as a common factor
[1, eq. (1.13), p. 7]:

$$
u(t) = K_p \left[ e(t) + \frac{1}{T_i} \int_0^t e(\tau) \, d\tau + T_d \frac{de(t)}{dt} \right]
$$

where $T_i = K_p / K_i$ is the integral time and $T_d = K_d / K_p$ is the
derivative time. The `isa_form` policy takes $K_p$, $T_i$ and $T_d$ in `kp`,
`ki` and `kd` and converts them to the parallel gains. This form makes the physical meaning of each parameter more
transparent: $T_i$ is the time for the integral term to repeat the proportional
action, and $T_d$ is the time by which the derivative term anticipates
the proportional action.

## Proportional Action

The proportional term produces output proportional to the current error
[1, Sec. 1.4.1, eq. (1.2), p. 3]:

$$
u_P(t) = K_p \, e(t)
$$

Increasing $K_p$ reduces steady-state error and increases the speed of
response, but excessive gain causes oscillation and eventually instability.
For most systems, proportional-only control leaves a residual steady-state
offset because the controller requires a non-zero error to produce a non-zero
output [3, Sec. 7.1, p. 338].

## Integral Action

The integral term accumulates past error, eliminating steady-state offset for
step disturbances and reference changes [1, Sec. 1.4.2, p. 5]. In the
Laplace domain, integral action contributes a pole at the origin, ensuring
infinite DC gain and hence zero steady-state error for constant inputs.

In discrete time with sampling period $T_s$, the backward-Euler approximation,
the library's default, is [1, Sec. 1.7, eq. (1.34), p. 14]:

$$
u_I[k] = u_I[k-1] + K_i \, T_s \, e[k]
$$

The trapezoidal (Tustin) approximation, selected by the `tustin` tag, provides better accuracy [5, Sec. 13.3,
eq. (13.9), p. 416]:

$$
u_I[k] = u_I[k-1] + \frac{K_i \, T_s}{2} \bigl(e[k] + e[k-1]\bigr)
$$

## Derivative Action

The derivative term anticipates future error based on its rate of change.
In practice, the derivative is rarely applied to the raw error: a setpoint step
would produce a derivative kick, and differentiation amplifies high-frequency
noise. The library differentiates $-y$ by default, and $c \cdot r - y$ with
`derivative_on_error` (see [Setpoint Weighting](#setpoint-weighting)). A
first-order low-pass filter with divisor $N$ gives the filtered derivative
[1, Sec. 1.6.1, eq. (1.20), p. 9]:

$$
D(s) = \frac{K_d \, s}{1 + s \, T_d / N}, \qquad T_d = K_d / K_p
$$

so the filter time constant is $T_f = T_d / N$, and the high-frequency gain is
limited to $K_p N$. The `deriv_filter` policy uses this default, or an explicit
$T_f$. The parameter $N$ determines the trade-off between derivative response speed
and noise attenuation. In most practical cases $N$ lies between 8 and 16
[1, Sec. 1.6.1, p. 9].

In discrete time the library uses the backward difference with period $T_s$,
here for the default derivative on the measurement:

$$
u_D[k] = \frac{T_f}{T_f + T_s} \, u_D[k-1] - \frac{K_d}{T_f + T_s} \bigl(y[k] - y[k-1]\bigr)
$$

## Anti-Windup

When the actuator saturates, the integral term continues accumulating error,
a phenomenon called integrator windup. This leads to large overshoots and
sluggish recovery when the setpoint changes [1, Sec. 3.2, pp. 35--36].

### Back-Calculation

Back-calculation feeds the difference between the limited and unlimited output
back into the integrator with tracking gain $K_b$
[1, Sec. 3.3.3, eq. (3.2), p. 39], here in the discrete form of [4, Algorithm 1]:

$$
u_I[k] = u_I[k-1] + K_i \, T_s \, e[k] + K_b \, T_s \bigl(u_{\text{sat}}[k-1] - u[k-1]\bigr)
$$

where $u$ is the command before rate limiting and saturation and
$u_{\text{sat}}$ the one applied.

The gain $K_b$ has units of inverse time. The default in this library is the
tracking-time-constant form $K_b = 1/T_t$ with $T_t = \sqrt{T_i T_d}$, which in
the internal parallel gains is $K_b = \sqrt{K_i/K_d}$. That rule cannot be used
without derivative action, so the library then takes the integral time,
$T_t = T_i$, giving the fallback $K_b = K_i/K_p = 1/T_i$
[1, Sec. 3.3.3, eqs. (3.3)--(3.4), p. 39].

### Clamping

Clamping stops integration when the output is saturated and the integrator
would make the saturation worse [1, Sec. 3.3.2, p. 38]. This is simpler than back-calculation but can be less smooth during
transitions.

In discrete time the saturation verdict depends on the output, and the output on
the integrator, so `anti_windup<clamping>` takes the verdict from the previous
cycle's constraint: the increment $\Delta u_I[k]$ is withheld when cycle $k-1$
was constrained and $\Delta u_I[k]$ has the sign of
$u[k-1] - u_{\text{sat}}[k-1]$. Taking the verdict from the previous sample's
saturation error is the standard scheme in reference implementations that form
the output from the pre-update integrator state [4, Algorithm 2], and it costs
at most one increment on the first constrained cycle. Algorithm 2 itself
withholds every increment after a constrained sample; the sign condition keeps
integrating when the integrator drives the output away from the limit, as the
last variant of [1, Sec. 3.3.2, p. 38] does. Comparing the increment's sign rather than the error's keeps the rule
correct for a reverse-acting controller with a negative $K_i$.

### Conditional Integration

`anti_windup<conditional_integration>` withholds the integral increment on
every cycle whose error magnitude exceeds a threshold, $|e[k]| > e_{\max}$,
whether or not the output saturates.

## Setpoint Weighting

Setpoint weighting separates the response to setpoint changes from the
response to disturbances. The proportional and derivative terms use weighted
setpoint signals [1, Secs. 1.6.2--1.6.3, eqs. (1.26) and (1.30), pp. 11--13]:

$$
e_P(t) = b \cdot r(t) - y(t), \qquad e_D(t) = c \cdot r(t) - y(t)
$$

where $b$ and $c$ are the proportional and derivative weights, usually in
$[0, 1]$; the library accepts any finite value. Setting $b < 1$ reduces
overshoot on setpoint steps without affecting disturbance rejection. Setting
$c = 0$ avoids derivative kicks from step changes in the reference. In the
library $c$ applies only with `derivative_on_error`; the default derivative
on $-y$ is the case $c = 0$. The integral acts on the unweighted error.

## Discrete-Time Implementation

A complete discrete PID implementation combines all the above elements.
The velocity (incremental) form computes the change in output rather than the
absolute output, avoiding the need to store the integral sum explicitly
[1, Sec. 1.7, eqs. (1.36)--(1.37), p. 14]:

$$
\Delta u[k] = u[k] - u[k-1]
$$

It eases bumpless switching between manual and automatic modes. In the
library's `velocity_form` the controller still keeps the running sum of its
increments to clamp it to the output limits, and a switch is bumpless when
that sum is re-seeded from the manual output through the tracking overload
of `compute`.

## Feed-Forward

Pure feedback control reacts only after an error appears. Feed-forward action
uses knowledge of the disturbance or setpoint trajectory to act preemptively
[1, Secs. 5.1--5.2, pp. 93--94]:

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

[4] M. Caparroz, K. Soltesz, T. Hagglund and J. L. Guzman, "Anti-Windup in PID
Control: Review, Analysis, and New Tuning Directions," arXiv:2606.01959, 2026.

[5] K. J. Astrom and T. Hagglund, "Advanced PID Control," ISA, 2006.
