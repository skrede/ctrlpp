# MPC Theory

Model Predictive Control (MPC) solves an optimal control problem at each time
step over a finite prediction horizon. It handles state and input constraints
naturally and provides a systematic way to trade off tracking performance
against control effort.

## Key Concepts

### Receding Horizon

At each time step k, MPC:

1. Measures (or estimates) the current state x[k]
2. Solves an optimisation problem over steps k to k+N
3. Applies only the first optimal input u*[k]
4. Advances to step k+1 and repeats

This "receding horizon" strategy re-plans at every step, providing implicit
robustness to model mismatch and disturbances.

### Cost Function

The standard quadratic cost penalises state deviation and control effort:

```
J = sum_{i=0}^{N-1} [ x[i]^T Q x[i] + u[i]^T R u[i] ] + x[N]^T Qf x[N]
```

- **Q** (state weight): penalises deviation from the target. Larger values
  mean tighter tracking.
- **R** (input weight): penalises control effort. Larger values mean smoother
  but slower control.
- **Qf** (terminal weight): penalises the final state. Choosing Qf as the
  solution to the DARE ensures closed-loop stability.

### Constraints

MPC naturally incorporates:

- **Input constraints**: actuator limits u_min <= u[i] <= u_max
- **State constraints**: safety limits x_min <= x[i] <= x_max
- **Terminal constraints**: x[N] in X_f (for stability guarantees)

Soft constraints with slack variables prevent infeasibility when hard state
constraints conflict with the dynamics.

### QP Formulation (Linear MPC)

For a linear system x[k+1] = A*x[k] + B*u[k], the MPC problem becomes a
quadratic program (QP). ctrlpp uses a sparse (non-condensed) formulation where
both states and inputs are decision variables:

```
minimise    (1/2) z^T H z + q^T z
subject to  A_eq z = b_eq       (dynamics constraints)
            l <= A_ineq z <= u   (state/input bounds)
```

The sparse formulation is more efficient for long horizons and enables
warm-starting between time steps.

### NLP Formulation (Nonlinear MPC)

For nonlinear dynamics x[k+1] = f(x[k], u[k]), the optimisation becomes a
nonlinear program (NLP). ctrlpp uses multiple shooting: each shooting segment
is a decision variable, and continuity constraints enforce that segments match
at boundaries. This provides better convergence properties than single
shooting.

### Terminal Constraints and Stability

Stability of MPC requires that the optimisation is recursively feasible and the
cost decreases over time. Sufficient conditions include:

- **Terminal cost Qf**: set to the solution of the discrete algebraic Riccati
  equation (DARE), making the terminal cost approximate the infinite-horizon
  cost
- **Terminal constraint set X_f**: an invariant set (e.g., ellipsoidal or
  polytopic) where the terminal LQR controller keeps the state feasible
  indefinitely

ctrlpp provides `terminal_ingredients()` to compute Qf and ellipsoidal
terminal sets from the system matrices.

## References

- **Rawlings, J. B., Mayne, D. Q., and Diehl, M.** *Model Predictive
  Control: Theory, Computation, and Design.* Nob Hill Publishing, 2nd ed.,
  2017. ISBN 978-0-9759377-3-0.
  The definitive MPC textbook. Covers stability, robustness, nonlinear MPC,
  and computational methods. Primary reference for ctrlpp's MPC
  implementation.

- **Mayne, D. Q., Rawlings, J. B., Rao, C. V., and Scokaert, P. O. M.**
  "Constrained Model Predictive Control: Stability and Optimality."
  *Automatica*, 36(6):789--814, 2000. DOI: 10.1016/S0005-1098(99)00214-9.
  Establishes the stability conditions (terminal cost + terminal set) used in
  ctrlpp's terminal_ingredients helper.

## Related API Pages

- [mpc](../mpc/mpc.md) -- linear MPC with sparse QP formulation
- [nmpc](../mpc/nmpc.md) -- nonlinear MPC with multiple shooting
- [osqp_solver](../mpc/osqp-solver.md) -- OSQP QP solver backend
- [nlopt_solver](../mpc/nlopt-solver.md) -- NLopt NLP solver backend
- [Solver Injection Guide](../guides/mpc/solver-injection.md) -- swapping
  solver backends
