# Nonlinear Model Predictive Control

Nonlinear Model Predictive Control (NMPC) extends the MPC framework to
systems with nonlinear dynamics and constraints. Instead of solving a
quadratic program at each time step, NMPC solves a nonlinear program (NLP),
enabling direct handling of nonlinear process models without linearization.

NMPC is essential for systems where linearization around an operating point
is insufficiently accurate, such as chemical reactors, aerospace vehicles,
and robotic manipulators operating over wide operating envelopes
[2, Sec. 8, pp. 391--417].

## Nonlinear MPC Formulation

For a discrete nonlinear system $x_{k+1} = f(x_k, u_k)$, the NMPC problem
at time $k$ is:

$$
\min_{x_{0:N}, u_{0:N-1}} \sum_{i=0}^{N-1} \ell(x_i, u_i) + V_f(x_N)
$$

$$
\text{subject to} \quad x_{i+1} = f(x_i, u_i), \quad i = 0, \ldots, N-1
$$

$$
g(x_i, u_i) \le 0, \quad i = 0, \ldots, N-1
$$

$$
h(x_N) \le 0
$$

$$
x_0 = \hat{x}_k
$$

where:

- $\ell(x, u)$ is the stage cost; ctrlpp's default is $\tfrac{1}{2}(x - x_r)^\top Q (x - x_r) + \tfrac{1}{2} u^\top R u$,
  with the input penalized about zero
- $V_f(x_N)$ is the terminal cost; ctrlpp's default is $\tfrac{1}{2}(x_N - x_r)^\top Q_f (x_N - x_r)$, with
  $Q_f = Q$ unless set
- $g(x, u) \le 0$ are path constraints (actuator limits, safety bounds)
- $h(x_N) \le 0$ are terminal constraints

ctrlpp's `nmpc` adds box bounds on the inputs and on the states of stages 1 to $N$, and an
optional symmetric bound on the input change per step.

## Direct Multiple Shooting

The most robust method for transcribing the NMPC problem into a finite-
dimensional NLP is direct multiple shooting [2, Sec. 3, pp. 395--400].

### Method

The prediction horizon is divided into $N$ intervals. On each interval $i$,
an independent initial value problem is solved:

$$
x_{i+1}^{\text{sim}} = F(x_i, u_i) = x_i + \int_{t_i}^{t_{i+1}} f(x(\tau), u_i) \, d\tau
$$

Both the states $x_0, \ldots, x_N$ and inputs $u_0, \ldots, u_{N-1}$ are
decision variables. Continuity between intervals is enforced by equality
constraints [2, Sec. 3.1, pp. 396--398]. ctrlpp takes $F$ to be the caller's
discrete-time dynamics callable itself, so any integration happens inside it:

$$
x_{i+1} - F(x_i, u_i) = 0, \quad i = 0, \ldots, N-1
$$

### Advantages

- **Robustness**: each shooting segment can be integrated independently,
  avoiding the numerical sensitivity of single shooting
- **Parallelism**: segment integrations are independent and could be run in
  parallel; ctrlpp evaluates them in sequence
- **Structure**: the NLP has a sparse, banded Jacobian that specialized
  solvers can exploit; ctrlpp's backends treat it as dense
- **Warm-starting**: the previous solution, shifted forward one stage, is the
  initial guess for the next time step, as ctrlpp's `nmpc` does

## NLP Solution Methods

### Sequential Quadratic Programming (SQP)

SQP solves the NLP by iterating over QP subproblems. At each iteration,
the cost is approximated by a quadratic model and the constraints by their
first-order Taylor expansions:

$$
\min_{\Delta z} \; \frac{1}{2} \Delta z^\top H_k \Delta z + \nabla f_k^\top \Delta z
$$

$$
\text{subject to} \quad \nabla c_k^\top \Delta z + c_k = 0
$$

where $H_k$ is an approximation to the Hessian of the Lagrangian. The QP
subproblem structure is identical to linear MPC, enabling reuse of efficient
QP solvers.

### Real-Time Iteration (RTI)

The real-time iteration scheme performs only a single SQP iteration per
control time step, using the warm-started solution from the previous step
[2, Sec. 5, pp. 405--410]. This trades optimality for computational speed:

1. **Preparation phase**: prepare the QP linearization using the predicted
   trajectory (can be done before the new measurement arrives)
2. **Feedback phase**: update only the initial state constraint with the
   new measurement and solve the QP

RTI converges over multiple control steps rather than within a single step,
achieving near-optimal performance with bounded computation time.

ctrlpp's `nmpc` does not implement RTI. Each control step runs the backend to
its own stopping criteria, from the shifted warm start; a small iteration cap
(`max_eval`) bounds the work per step but does not split it into preparation
and feedback phases.

## Constraint Handling

### Path Constraints

General nonlinear path constraints $g(x_i, u_i) \le 0$ enter the NLP
directly. The NLP solver handles them via active-set or interior-point methods.
ctrlpp ships no interior-point backend: argmin's SQP policies solve active-set
QP subproblems, and the other backends handle the rows by moving asymptotes,
derivative-free models or an augmented Lagrangian.

### Soft Constraints

As in linear MPC, infeasibility can occur when hard constraints conflict with
the dynamics. Soft constraints with L1 penalty relaxation prevent solver
failure [1, Sec. 12.6, pp. 264--266]:

$$
g(x_i, u_i) - \epsilon_i \le 0, \quad \epsilon_i \ge 0
$$

with additional cost $\rho \sum_i \epsilon_i$.

ctrlpp softens the path and terminal constraints this way, with one
non-negative slack per constraint and stage and a per-constraint weight
$\rho$ (`path_penalty`, `terminal_penalty`); the state and input bounds stay
hard. This differs from ctrlpp's linear `mpc`, which softens its state bounds
with a quadratic penalty instead.

## Stability

Stability analysis for NMPC uses the same sufficient conditions as linear MPC;
for a nonlinear model they are met locally, from the linearization at the
origin [3, Secs. 3.2--3.3, pp. 796--797; Sec. 3.7.4.3, p. 801]:

1. **Terminal cost**: $V_f(x)$ must be a local control Lyapunov function
   in a neighborhood of the origin
2. **Terminal constraint**: $x_N$ must be constrained to a region where
   a local stabilizing controller exists
3. **Continuity**: the value function must be continuous, the standing
   assumption of the analysis

When these conditions hold, the NMPC value function is a Lyapunov function
and the closed loop is asymptotically stable. ctrlpp computes no terminal
ingredients for `nmpc`: the caller supplies the terminal cost (`Qf` or
`terminal_cost_cb`) and the terminal constraint (`terminal_constraint_cb`).

## Comparison with Linear MPC

| Aspect           | Linear MPC           | NMPC                        |
| ---------------- | -------------------- | --------------------------- |
| Plant model      | Linear               | Nonlinear                   |
| Optimization     | QP (convex)          | NLP (non-convex)            |
| Global optimum   | Guaranteed           | Local only                  |
| Computation      | Milliseconds         | Milliseconds to seconds     |
| Warm-starting    | Very effective       | Critical for convergence    |
| Stability proof  | Well-established     | Requires terminal ingredients|

## References

[1] F. Borrelli, A. Bemporad, and M. Morari, "Predictive Control for Linear
and Hybrid Systems," Cambridge University Press, 2017.

[2] M. Diehl, H. J. Ferreau, and N. Haverbeke, "Efficient Numerical Methods
for Nonlinear MPC and Moving Horizon Estimation," in Nonlinear Model
Predictive Control, LNCIS vol. 384, Springer, 2009, pp. 391--417.

[3] D. Q. Mayne, J. B. Rawlings, C. V. Rao, and P. O. M. Scokaert,
"Constrained Model Predictive Control: Stability and Optimality," Automatica,
vol. 36, no. 6, pp. 789--814, 2000.

[4] J. B. Rawlings, D. Q. Mayne, and M. Diehl, "Model Predictive Control:
Theory, Computation, and Design," 2nd ed., Nob Hill Publishing, 2017.
