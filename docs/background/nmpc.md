# Nonlinear Model Predictive Control

Nonlinear Model Predictive Control (NMPC) extends the MPC framework to
systems with nonlinear dynamics and constraints. Instead of solving a
quadratic program at each time step, NMPC solves a nonlinear program (NLP),
enabling direct handling of nonlinear process models without linearisation
[1, Ch. 12, pp. 221--250].

NMPC is essential for systems where linearisation around an operating point
is insufficiently accurate, such as chemical reactors, aerospace vehicles,
and robotic manipulators operating over wide operating envelopes
[2, Sec. 8, pp. 391--417].

## Nonlinear MPC Formulation

For a discrete nonlinear system $x_{k+1} = f(x_k, u_k)$, the NMPC problem
at time $k$ is [1, Sec. 12.1, pp. 222--226]:

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

- $\ell(x, u)$ is the stage cost (often quadratic: $(x - x_r)^\top Q (x - x_r) + u^\top R u$)
- $V_f(x_N)$ is the terminal cost
- $g(x, u) \le 0$ are path constraints (actuator limits, safety bounds)
- $h(x_N) \le 0$ are terminal constraints

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
constraints [2, Sec. 3.1, pp. 396--398]:

$$
x_{i+1} - F(x_i, u_i) = 0, \quad i = 0, \ldots, N-1
$$

### Advantages

- **Robustness**: each shooting segment can be integrated independently,
  avoiding the numerical sensitivity of single shooting
- **Parallelism**: segment integrations are independent and can be
  parallelised
- **Structure**: the NLP has a sparse, banded Jacobian that specialised
  solvers can exploit
- **Warm-starting**: the previous solution provides an excellent initial
  guess for the next time step by shifting the trajectory forward

## NLP Solution Methods

### Sequential Quadratic Programming (SQP)

SQP solves the NLP by iterating over QP subproblems. At each iteration,
the cost is approximated by a quadratic model and the constraints by their
first-order Taylor expansions [1, Sec. 12.3, pp. 232--236]:

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

1. **Preparation phase**: prepare the QP linearisation using the predicted
   trajectory (can be done before the new measurement arrives)
2. **Feedback phase**: update only the initial state constraint with the
   new measurement and solve the QP

RTI converges over multiple control steps rather than within a single step,
achieving near-optimal performance with bounded computation time.

## Constraint Handling

### Path Constraints

General nonlinear path constraints $g(x_i, u_i) \le 0$ enter the NLP
directly. The NLP solver handles them via active-set or interior-point methods
[1, Sec. 12.4, pp. 237--240].

### Soft Constraints

As in linear MPC, infeasibility can occur when hard state constraints conflict
with the dynamics. Soft constraints with L1 penalty relaxation prevent
solver failure [1, Sec. 12.4, pp. 240--242]:

$$
g(x_i, u_i) - \epsilon_i \le 0, \quad \epsilon_i \ge 0
$$

with additional cost $\rho \sum_i \epsilon_i$.

## Stability

Stability analysis for NMPC follows the same framework as linear MPC but with
stronger assumptions [3, Sec. 5, pp. 800--808]:

1. **Terminal cost**: $V_f(x)$ must be a local control Lyapunov function
   in a neighbourhood of the origin
2. **Terminal constraint**: $x_N$ must be constrained to a region where
   a local stabilising controller exists
3. **Continuity**: the stage cost and dynamics must be continuous

When these conditions hold, the NMPC value function is a Lyapunov function
and the closed loop is asymptotically stable.

## Comparison with Linear MPC

| Aspect           | Linear MPC           | NMPC                        |
| ---------------- | -------------------- | --------------------------- |
| Plant model      | Linear               | Nonlinear                   |
| Optimisation     | QP (convex)          | NLP (non-convex)            |
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
