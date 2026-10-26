# Model Predictive Control

Model Predictive Control (MPC) is an optimization-based control strategy that
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
2. Solve an optimization problem over a horizon of $N$ steps
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

- $Q \succeq 0$ penalizes state deviation from the target
- $R \succ 0$ penalizes control effort
- $Q_f \succeq 0$ penalizes the terminal state
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
block-sparse, enabling efficient specialized solvers. The sparse formulation
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

A terminal constraint $x_N \in \mathcal{X}_f$ asks the horizon to end in a set
where the LQR law $u = K x$, with $K = -(R + B^\top P B)^{-1} B^\top P A$ and
$P$ the DARE solution, keeps the state and the input admissible from then on
[3, Sec. 4, pp. 800--808]. The library offers three sets, and only the first is
positively invariant under that law once it is encoded in the QP:

- **Maximal output admissible set**
  $\mathcal{O}_\infty = \{x : A_K^k x \in \mathcal{X},\ K A_K^k x \in \mathcal{U}\ \text{for all}\ k \ge 0\}$
  with $A_K = A + B K$: the largest set inside $\{x \in \mathcal{X},\ K x \in \mathcal{U}\}$
  that $A_K$ maps into itself [4]. When $A_K$ is strictly stable, the finitely
  bounded outputs observe $A_K$, and the bounds form a compact polytope with the
  origin in its interior, $\mathcal{O}_\infty$ is a compact polytope described by
  finitely many of the rows $C A_K^k$, $C = [I;\ K]$ [4], [5].
  `terminal_ingredients` returns it as `invariant_set`, and
  `compute_output_admissible_set` builds it for any stable closed loop with
  $n_x \le 4$, by vertex enumeration and with a step bound from a Lyapunov
  contraction of $A_K$ [5, Theorem 5].
- **Ellipsoid** $\{x : x^\top P x \le \alpha\}$, with $\alpha$ the largest level
  inside the input and state faces (`terminal_ingredients`' `set`). A QP cannot
  hold a quadric, so the controller encodes the ellipsoid's inscribed box in
  $P$'s eigenbasis. That box is not positively invariant under $u = K x$ and
  carries no recursive-feasibility guarantee.
- **Robust polytope** (`compute_polytopic_invariant_set`): the states that stay
  inside the state box under every input in $\mathcal{U}$, a robust
  control-invariant set. It is a valid object for robust MPC, but it is not the
  terminal set of the LQR law.

With the DARE terminal cost, the same gain $K$, $\mathcal{X}_f = \mathcal{O}_\infty$,
a nominal model, regulation to the origin, and hard state constraints, the
closed loop is recursively feasible and the origin is asymptotically stable
[3]. The controller softens state constraints by default
(`hard_state_constraints = false`). With softened state constraints the
optimization admits state-bound violations at a penalty, so the feasibility and
cost-decrease argument changes and the statement does not apply as written.
Tracking a nonzero reference is also outside it.

## Reference Tracking

For setpoint tracking, the cost function penalizes deviation from a target
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

[4] E. G. Gilbert and K. T. Tan, "Linear Systems with State and Control
Constraints: The Theory and Application of Maximal Output Admissible Sets,"
IEEE Transactions on Automatic Control, vol. 36, no. 9, pp. 1008--1020, 1991.

[5] H. R. Ossareh and I. Kolmanovsky, "On Complexity Bounds for the Maximal
Admissible Set of Linear Time-Invariant Systems," arXiv:2302.02246, 2023.
