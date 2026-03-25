# MPC Theory

Model Predictive Control (MPC) solves an optimal control problem at each time
step over a finite prediction horizon. It handles state and input constraints
naturally and provides a systematic way to trade off tracking performance
against control effort.

## Receding Horizon

At each time step $k$, MPC:

1. Measures (or estimates) the current state $x_k$
2. Solves an optimisation problem over steps $k$ to $k + N$
3. Applies only the first optimal input $u_0^\star$
4. Advances to step $k + 1$ and repeats

This "receding horizon" strategy re-plans at every step, providing implicit
robustness to model mismatch and disturbances.

## Optimisation Problem

### Linear MPC

For a linear time-invariant system $x_{k+1} = A x_k + B u_k$, the MPC
problem is:

$$
\min_{u_0, \ldots, u_{N-1}} \sum_{i=0}^{N-1} \left( x_i^\top Q \, x_i + u_i^\top R \, u_i \right) + x_N^\top Q_f \, x_N
$$

$$
\text{subject to} \quad x_{i+1} = A \, x_i + B \, u_i, \quad i = 0, \ldots, N-1
$$

$$
u_{\min} \le u_i \le u_{\max}, \quad x_{\min} \le x_i \le x_{\max}
$$

$$
x_N \in \mathcal{X}_f
$$

where:

- $Q \succeq 0$ (state weight) penalises deviation from the target
- $R \succ 0$ (input weight) penalises control effort
- $Q_f \succeq 0$ (terminal weight) penalises the final state
- $\mathcal{X}_f$ is the terminal constraint set

### QP Standard Form

For the linear case, the optimisation reduces to a quadratic program (QP).
ctrlpp uses a sparse (non-condensed) formulation where both states and inputs
are decision variables. Stacking $z = [x_0^\top, u_0^\top, x_1^\top, u_1^\top, \ldots, x_N^\top]^\top$:

$$
\min_z \; \frac{1}{2} z^\top H z + q^\top z
$$

$$
\text{subject to} \quad A_{\text{eq}} \, z = b_{\text{eq}} \quad \text{(dynamics)}
$$

$$
l \le A_{\text{ineq}} \, z \le u \quad \text{(bounds)}
$$

The sparse formulation preserves problem structure, is more efficient for long
horizons, and enables warm-starting between time steps.

### Nonlinear MPC

For nonlinear dynamics $x_{k+1} = f(x_k, u_k)$, the optimisation becomes a
nonlinear program (NLP):

$$
\min_{x_{0:N}, u_{0:N-1}} \sum_{i=0}^{N-1} \ell(x_i, u_i) + V_f(x_N)
$$

$$
\text{subject to} \quad x_{i+1} = f(x_i, u_i)
$$

$$
g(x_i, u_i) \le 0, \quad h(x_N) \le 0
$$

ctrlpp uses multiple shooting: each shooting segment is a decision variable,
and continuity constraints enforce that segments match at boundaries. This
provides better convergence properties than single shooting.

## Constraints

MPC naturally incorporates:

- **Input constraints**: actuator limits $u_{\min} \le u_i \le u_{\max}$
- **State constraints**: safety limits $x_{\min} \le x_i \le x_{\max}$
- **Path constraints**: general nonlinear $g(x_i, u_i) \le 0$
- **Terminal constraints**: $x_N \in \mathcal{X}_f$ (for stability)

Soft constraints with slack variables prevent infeasibility when hard state
constraints conflict with the dynamics. ctrlpp uses L1 soft constraints by
default with configurable penalty weight.

## Terminal Constraints and Stability

Stability of MPC requires that the optimisation is recursively feasible and
the cost decreases over time. Sufficient conditions (Mayne et al., 2000):

1. **Terminal cost** $Q_f$: set to the solution of the discrete algebraic
   Riccati equation (DARE):

$$
Q_f = A^\top Q_f A - A^\top Q_f B (R + B^\top Q_f B)^{-1} B^\top Q_f A + Q
$$

   This makes the terminal cost approximate the infinite-horizon LQR cost.

2. **Terminal constraint set** $\mathcal{X}_f$: an invariant set where the
   LQR controller $u = -K_\infty x$ keeps the state feasible. ctrlpp
   supports ellipsoidal sets $\{x : x^\top P^{-1} x \le 1\}$ and polytopic
   sets.

ctrlpp provides `terminal_ingredients()` to compute $Q_f$, $K_\infty$, and
the ellipsoidal terminal set from the system matrices.

## References

- Rawlings, J. B., Mayne, D. Q., and Diehl, M. (2017). *Model Predictive
  Control: Theory, Computation, and Design.* Nob Hill Publishing, 2nd ed.
  [`rawlings2017`]
  The definitive MPC textbook. Covers stability, robustness, nonlinear MPC,
  and computational methods. Primary reference for ctrlpp's MPC
  implementation.

- Mayne, D. Q., Rawlings, J. B., Rao, C. V., and Scokaert, P. O. M. (2000).
  "Constrained Model Predictive Control: Stability and Optimality."
  *Automatica*, 36(6):789--814. [`mayne2000`]
  Establishes the stability conditions (terminal cost + terminal set) used in
  ctrlpp's `terminal_ingredients` helper.

## Related API Pages

- [mpc](../mpc/mpc.md) -- linear MPC with sparse QP formulation
- [nmpc](../mpc/nmpc.md) -- nonlinear MPC with multiple shooting
- [osqp_solver](../mpc/osqp-solver.md) -- OSQP QP solver backend
- [nlopt_solver](../mpc/nlopt-solver.md) -- NLopt NLP solver backend
- [Solver Injection Guide](../guides/mpc/solver-injection.md) -- swapping
  solver backends
