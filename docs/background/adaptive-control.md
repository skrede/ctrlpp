# Adaptive Control

Adaptive control addresses the fundamental challenge of controlling systems
whose parameters are unknown or change over time. Unlike fixed-gain controllers
such as PID or LQR that require accurate plant models, adaptive controllers
adjust their parameters online to maintain desired performance despite model
uncertainty. The key insight is that while the exact plant parameters may be
unknown, the plant structure (order, relative degree, sign of the input gain)
is assumed known [1, Ch. 8, pp. 309--312].

The importance of robustification in adaptive control was dramatically
illustrated by the Rohrs counterexample (1985), which showed that even small
unmodeled dynamics or bounded disturbances could drive an unrobustified
adaptive controller unstable. This result motivated the development of
robust adaptive laws including dead-zone modification, sigma-modification,
and e-modification [1, Sec. 8.6, pp. 367--380].

## MRAC

Model Reference Adaptive Control (MRAC) is the most widely studied direct
adaptive control approach. The controller adjusts its parameters so that the
closed-loop plant behavior matches a designer-specified reference model,
without explicitly identifying the plant parameters.

### Problem Statement

Consider an unknown discrete-time SISO plant with known structure:

$$
x_{k+1} = a_p \, x_k + b_p \, u_k
$$

where the coefficients $a_p$ and $b_p$ are unknown but the sign of $b_p$ is
assumed known. The designer specifies a stable reference model that captures
the desired closed-loop behavior:

$$
x_{m,k+1} = a_m \, x_{m,k} + b_m \, r_k
$$

where $a_m$ and $b_m$ are chosen such that $|a_m| < 1$ (stable pole inside the
unit circle) and $r_k$ is the external reference signal. The control objective
is to find a control law $u_k$ such that the tracking error:

$$
e_k = x_k - x_{m,k}
$$

converges to zero as $k \to \infty$, despite the unknown plant parameters.

### Adaptation Law

The control law takes the form of a linear state-feedback plus reference
feedforward with adaptive gains [1, Sec. 8.4, pp. 339--346]:

$$
u_k = \theta_{x,k} \, x_k + \theta_{r,k} \, r_k
$$

where $\theta_{x,k}$ and $\theta_{r,k}$ are the adaptive parameters. If the
plant were known, the ideal parameters would be
$\theta_x^* = (a_m - a_p) / b_p$ and $\theta_r^* = b_m / b_p$, which make
the closed-loop dynamics identical to the reference model. Since $a_p$ and
$b_p$ are unknown, the parameters must be adapted online.

The Lyapunov-based adaptation law is derived by constructing a positive
definite function of the tracking error and parameter errors, then choosing
update rules that guarantee the function is non-increasing
[1, Sec. 8.4, pp. 342--345]. The resulting discrete-time update equations are:

$$
\theta_{x,k+1} = \theta_{x,k} - \gamma \, \mathrm{sgn}(b_p) \, e_k \, x_k
$$

$$
\theta_{r,k+1} = \theta_{r,k} - \gamma \, \mathrm{sgn}(b_p) \, e_k \, r_k
$$

where $\gamma > 0$ is the adaptation gain. Larger $\gamma$ produces faster
adaptation but can cause oscillation; smaller $\gamma$ gives smoother
convergence but slower tracking. The sign of $b_p$ ensures the adaptation
moves in the correct direction regardless of the unknown input gain magnitude.

The Lyapunov argument guarantees that $e_k \to 0$ and the parameters remain
bounded, but does not guarantee convergence of $\theta_x$ and $\theta_r$ to
their ideal values unless a persistence of excitation condition is satisfied
[2, Ch. 4, pp. 153--162].

### Robustification

The Rohrs counterexample (1985) demonstrated that the basic MRAC law above
can become unstable in the presence of unmodeled dynamics, bounded
disturbances, or measurement noise. Even small perturbations not captured by
the assumed plant structure can cause unbounded parameter drift, leading to
instability. This motivated the development of several robustification
strategies that trade a small amount of tracking accuracy for guaranteed
parameter boundedness.

#### Dead-Zone Modification

Dead-zone modification freezes adaptation when the tracking error is small,
preventing parameter drift driven by noise or modeling errors rather than
genuine tracking deviation [1, Sec. 8.6.1, pp. 368--371]:

$$
\theta_{x,k+1} = \begin{cases}
\theta_{x,k} - \gamma \, \mathrm{sgn}(b_p) \, e_k \, x_k & \text{if } |e_k| > \epsilon_0 \\
\theta_{x,k} & \text{otherwise}
\end{cases}
$$

$$
\theta_{r,k+1} = \begin{cases}
\theta_{r,k} - \gamma \, \mathrm{sgn}(b_p) \, e_k \, r_k & \text{if } |e_k| > \epsilon_0 \\
\theta_{r,k} & \text{otherwise}
\end{cases}
$$

The threshold $\epsilon_0$ must be set larger than the worst-case error due to
unmodeled effects. The tracking error converges to a ball of radius
$\epsilon_0$ rather than to zero.

#### Sigma-Modification

Sigma-modification adds a leaky integrator term that continuously pulls
the adaptive parameters toward zero [3, Ch. 8, pp. 412--418]:

$$
\theta_{x,k+1} = \theta_{x,k} - \gamma \left( \mathrm{sgn}(b_p) \, e_k \, x_k + \sigma \, \theta_{x,k} \right)
$$

$$
\theta_{r,k+1} = \theta_{r,k} - \gamma \left( \mathrm{sgn}(b_p) \, e_k \, r_k + \sigma \, \theta_{r,k} \right)
$$

where $\sigma > 0$ is the leakage coefficient. This guarantees bounded
parameters regardless of disturbances. The trade-off is a steady-state bias:
the parameters are pulled away from their ideal values, introducing a
non-zero steady-state tracking error proportional to $\sigma$.

#### e-Modification

e-modification scales the leakage term by the magnitude of the tracking
error, so that leakage is active only when the error is significant
[2, Ch. 8, pp. 385--392]:

$$
\theta_{x,k+1} = \theta_{x,k} - \gamma \left( \mathrm{sgn}(b_p) \, e_k \, x_k + \delta \, |e_k| \, \theta_{x,k} \right)
$$

$$
\theta_{r,k+1} = \theta_{r,k} - \gamma \left( \mathrm{sgn}(b_p) \, e_k \, r_k + \delta \, |e_k| \, \theta_{r,k} \right)
$$

where $\delta > 0$ is the e-modification gain. When the tracking error is
small (indicating good adaptation), the leakage is negligible and the
parameters settle near their ideal values. When the error is large (indicating
possible instability), the leakage becomes significant and pulls parameters
back. This makes e-modification preferred when the tracking error is a
reliable indicator of adaptation quality.

### ctrlpp Implementation

The library provides `mrac_controller<Scalar, NX, NU, Robustification>` as
a stateful controller with compile-time robustification policy selection.
The robustification mode is a template parameter (one of
`no_robustification`, `dead_zone`, `sigma_modification`, or
`e_modification`), enabling zero-overhead dispatch via `if constexpr`.
The `evaluate(x, r)` method computes the control signal and updates the
adaptive parameters in a single call, matching the one-step-per-call
convention used throughout ctrlpp.

## L1 Adaptive Control

L1 adaptive control is a more recent architecture that decouples the
adaptation rate from the robustness guarantees, addressing a fundamental
limitation of classical MRAC where fast adaptation can excite unmodeled
dynamics [4, Ch. 1, pp. 1--12].

The L1 architecture consists of three components:

- **State predictor**: a model that predicts the plant state based on current
  parameter estimates, producing a prediction error that drives adaptation.
- **Adaptation law**: a fast (potentially arbitrarily fast) parameter update
  that minimizes the prediction error. Unlike MRAC, the adaptation rate does
  not directly affect the closed-loop bandwidth.
- **Low-pass filter**: a control filter $C(z)$ that limits the bandwidth of
  the control signal, preventing high-frequency content from exciting
  unmodeled dynamics. The filter bandwidth is the sole robustness knob,
  separating it cleanly from the adaptation speed.

This architecture guarantees bounded transient performance (the tracking
error stays within a computable bound from the first time step) and allows
systematic trade-off between performance and robustness through the filter
design [4, Ch. 2, pp. 13--42].

Implementation details for L1 adaptive control are planned for Phase 34.

## References

[1] J.-J. E. Slotine and W. Li, "Applied Nonlinear Control," Prentice Hall,
1991.

[2] K. S. Narendra and A. M. Annaswamy, "Stable Adaptive Systems," Dover,
2005.

[3] P. A. Ioannou and J. Sun, "Robust Adaptive Control," Dover, 2012.

[4] N. Hovakimyan and C. Cao, "L1 Adaptive Control Theory: Guaranteed
Robustness with Fast Adaptation," SIAM, 2010.
