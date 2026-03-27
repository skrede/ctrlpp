# Model Predictive Control

Model Predictive Control (MPC) is an optimisation-based control strategy that
solves a finite-horizon optimal control problem at each time step, applies
only the first element of the optimal input sequence, and repeats at the next
step. This receding horizon approach handles state and input constraints
naturally and provides a systematic way to balance tracking performance
against control effort [1, Ch. 11, pp. 187--220].

MPC has become the dominant advanced control strategy in the process
industries and is increasingly used in aerospace, automotive, and robotics
applications where constraint satisfaction is critical
[2, Preface, pp. ii--iii].

## Receding Horizon Principle

At each time step $k$, MPC performs four operations
[1, Sec. 11.1, pp. 188--190]:

1. Measure (or estimate) the current state $x_k$
2. Solve an optimisation problem over a horizon of $N$ steps
3. Apply only the first optimal input $u_0^\star$ to the plant
4. Advance to step $k + 1$ and repeat

This "receding horizon" strategy re-plans at every step, providing implicit
feedback and robustness to model mismatch and disturbances. The horizon
slides forward in time, always looking $N$ steps ahead from the current
state.

## Linear MPC Formulation

For a discrete linear time-invariant system $x_{k+1} = A x_k + B u_k$ with
output $y_k = C x_k$, the MPC problem at time $k$ is
[1, Sec. 11.2, pp. 191--198]:

$$
\min_{u_0, \ldots, u_{N-1}} \sum_{i=0}^{N-1} \left( x_i^\top Q \, x_i + u_i^\top R \, u_i \right) + x_N^\top Q_f \, x_N
$$

$$
\text{subject to} \quad x_{i+1} = A \, x_i + B \, u_i, \quad i = 0, \ldots, N-1
$$

$$
u_{\min} \le u_i \le u_{\max}, \quad x_{\min} \le x_i \le x_{\max}
$$

where:

- $Q \succeq 0$ penalises state deviation from the target
- $R \succ 0$ penalises control effort
- $Q_f \succeq 0$ penalises the terminal state
- The inequality constraints encode actuator limits and safety bounds

## QP Standard Form

For linear systems with quadratic cost, the MPC problem reduces to a
quadratic program (QP). Two formulations exist [1, Sec. 11.3, pp. 199--205]:

### Dense (Condensed) Formulation

Eliminate the state variables using the dynamics equations, expressing all
states as functions of $u_0, \ldots, u_{N-1}$ and the initial state $x_0$.
The decision variable is $U = [u_0^\top, \ldots, u_{N-1}^\top]^\top$:

$$
\min_U \; \frac{1}{2} U^\top H \, U + f^\top U
$$

This has $N \cdot n_u$ decision variables. Efficient for short horizons
and small input dimensions.

### Sparse (Non-Condensed) Formulation

Keep both states and inputs as decision variables. Stacking
$z = [x_0^\top, u_0^\top, x_1^\top, u_1^\top, \ldots, x_N^\top]^\top$:

$$
\min_z \; \frac{1}{2} z^\top H \, z + q^\top z
$$

$$
\text{subject to} \quad A_{\text{eq}} \, z = b_{\text{eq}} \quad \text{(dynamics)}
$$

$$
l \le z \le u \quad \text{(bounds)}
$$

This has $N \cdot (n_x + n_u) + n_x$ decision variables but the matrices are
block-sparse, enabling efficient specialised solvers. The sparse formulation
preserves problem structure, is more efficient for long horizons, and enables
warm-starting between time steps [1, Sec. 11.3, p. 203].

## Constraints

MPC naturally incorporates [1, Sec. 11.4, pp. 206--210]:

- **Input constraints**: actuator limits $u_{\min} \le u_i \le u_{\max}$
- **State constraints**: safety limits $x_{\min} \le x_i \le x_{\max}$
- **Rate constraints**: $\Delta u_{\min} \le u_i - u_{i-1} \le \Delta u_{\max}$

### Soft Constraints

Hard state constraints can render the QP infeasible if the dynamics cannot
satisfy them from the current state. Soft constraints add slack variables
with penalty weights to maintain feasibility
[1, Sec. 11.4, pp. 210--212]:

$$
x_{\min} - \epsilon_i \le x_i \le x_{\max} + \epsilon_i, \qquad \epsilon_i \ge 0
$$

with an additional cost term $\rho \sum_i \lVert \epsilon_i \rVert_1$ where
$\rho$ is the penalty weight. L1 penalties produce exact constraint
satisfaction when the penalty is sufficiently large.

## Terminal Cost and Stability

Without a terminal cost, MPC with finite horizon provides no stability
guarantees. The seminal work of Mayne et al. [3, pp. 789--814] established
sufficient conditions for stability:

### Terminal Cost $Q_f$

Setting $Q_f$ to the solution of the discrete algebraic Riccati equation
(DARE) makes the terminal cost approximate the infinite-horizon LQR cost
[3, Sec. 3, pp. 795--800]:

$$
Q_f = A^\top Q_f A - A^\top Q_f B (R + B^\top Q_f B)^{-1} B^\top Q_f A + Q
$$

### Terminal Constraint Set $\mathcal{X}_f$

A terminal constraint $x_N \in \mathcal{X}_f$ ensures the state reaches an
invariant region where the unconstrained LQR controller $u = -K_\infty x$
keeps the state feasible [3, Sec. 4, pp. 800--808]:

- **Ellipsoidal sets**: $\{x : x^\top P^{-1} x \le 1\}$ where $P$ is the
  DARE solution
- **Polytopic sets**: $\{x : F x \le g\}$ computed from the maximal positive
  invariant set under $u = -K_\infty x$

Together, the DARE terminal cost and a positive invariant terminal set
guarantee recursive feasibility and asymptotic stability of the origin.

## Reference Tracking

For setpoint tracking, the cost function penalises deviation from a target
$(x_{\text{ref}}, u_{\text{ref}})$ rather than from the origin
[1, Sec. 11.5, pp. 213--215]:

$$
\min \sum_{i=0}^{N-1} \left( (x_i - x_{\text{ref}})^\top Q (x_i - x_{\text{ref}}) + (u_i - u_{\text{ref}})^\top R (u_i - u_{\text{ref}}) \right) + \ldots
$$

The reference pair must satisfy the steady-state equation
$x_{\text{ref}} = A x_{\text{ref}} + B u_{\text{ref}}$.

## References

[1] F. Borrelli, A. Bemporad, and M. Morari, "Predictive Control for Linear
and Hybrid Systems," Cambridge University Press, 2017.

[2] J. B. Rawlings, D. Q. Mayne, and M. Diehl, "Model Predictive Control:
Theory, Computation, and Design," 2nd ed., Nob Hill Publishing, 2017.

[3] D. Q. Mayne, J. B. Rawlings, C. V. Rao, and P. O. M. Scokaert,
"Constrained Model Predictive Control: Stability and Optimality," Automatica,
vol. 36, no. 6, pp. 789--814, 2000.
