# Nonlinear Moving Horizon Estimation

Nonlinear Moving Horizon Estimation (NMHE) extends the MHE framework to
systems with nonlinear dynamics and measurement models. Instead of solving a
quadratic program, NMHE solves a nonlinear program (NLP) over a sliding
window of recent measurements, enabling state estimation for general nonlinear
constrained systems [1, Sec. 1, pp. 246--248].

NMHE is the estimation dual of nonlinear MPC: both solve NLPs at each time
step and both use receding horizons. In ctrlpp, `nmhe` and the NMPC
controllers take the same NLP backends through the `nlp_solver` concept.

## NMHE Formulation

For a discrete nonlinear system $x_{i+1} = f(x_i, u_i) + w_i$ with
measurements $z_i = h(x_i) + v_i$, the NMHE problem at time $k$ is
[1, Sec. 1, p. 248]:

$$
\min_{x_{k-N}, \, w_{k-N:k-1}} \;
\lVert x_{k-N} - \bar{x}_{k-N} \rVert^2_{P_{\text{arr}}^{-1}}
+ \sum_{i=k-N}^{k-1} \lVert w_i \rVert^2_{Q^{-1}}
+ \sum_{i=k-N}^{k} \lVert z_i - h(x_i) \rVert^2_{R^{-1}}
$$

$$
\text{subject to} \quad x_{i+1} = f(x_i, u_i) + w_i
$$

$$
x_{\min} \le x_i \le x_{\max}
$$

The key differences from linear MHE are:
- The dynamics constraint uses the nonlinear model $f$ instead of linear $Ax + Bu$
- The measurement residual uses the nonlinear observation $h(x)$ instead of $Hx$
- The resulting optimization is a nonlinear program, not a QP

ctrlpp's `nmhe` uses the same convention as its `mhe`: every state is a
decision variable, the process noise is the residual $x_{i+1} - f(x_i, u_i)$
inside the cost, every term carries a factor $\tfrac{1}{2}$, and the arrival
term is scaled by `arrival_cost_weight`. Its box bounds are hard. It adds a
hard bound on the measurement residuals and a per-stage path constraint
$g(x_i) \le 0$ that can be softened (see [State Constraints](#state-constraints)).

## NLP Structure

### Decision Variables

The NLP decision variables are the initial state $x_{k-N}$ and the process
noise sequence $w_{k-N}, \ldots, w_{k-1}$. The remaining states are
determined implicitly by the dynamics constraints
[1, Sec. 1, p. 248]:

$$
x_{i+1} = f(x_i, u_i) + w_i, \quad i = k-N, \ldots, k-1
$$

Alternatively, all states can be treated as explicit decision variables
(multiple shooting formulation) with continuity constraints, providing
better convergence properties [2, Sec. 3, pp. 395--400].

ctrlpp's `nmhe` keeps all states as explicit decision variables but has no
continuity constraints: a continuity equality would pin every $w_i$ to zero,
so the dynamics enter only through the weighted residual in the cost.

### Gradient and Hessian

The NLP gradient requires the Jacobians of $f$ and $h$ with respect to the
state. These can be computed analytically or by finite differences. The Hessian
of the Lagrangian can be approximated using a Gauss-Newton method, which
exploits the least-squares structure of the cost function
[2, Sec. 4, pp. 400--405].

ctrlpp's `nmhe` uses neither analytic Jacobians nor a Gauss-Newton Hessian.
It computes the cost gradient by central differences of the whole cost,
leaves the constraint Jacobian to the backend's finite differences, and
relies on the backend's own Hessian approximation.

## Arrival Cost for Nonlinear Systems

The arrival cost in NMHE must summarize information from all measurements
before the estimation window. For nonlinear systems, the exact arrival cost
is intractable, so approximations are used [1, Sec. 1, p. 248]:

### EKF-Based Arrival Cost

A companion EKF runs alongside the NMHE. At each window shift, the EKF
provides:
- The prior mean $\bar{x}_{k-N}$ from its state estimate
- The arrival cost weighting $P_{\text{arr}}$ from its error covariance

This is the most common approach and provides adequate approximations for
mildly nonlinear systems [2, Sec. 4.3, pp. 403--405].

### Particle-Based Arrival Cost

For strongly nonlinear or multi-modal systems, a particle filter
approximation of the arrival cost can be used. This captures non-Gaussian
prior information but is computationally more expensive. ctrlpp implements
only the EKF-based arrival cost.

## State Constraints

NMHE handles state constraints directly in the NLP formulation. Common
constraint types [1, Sec. 1, p. 246]:

- **Box constraints**: $x_{\min} \le x_i \le x_{\max}$ (non-negative
  concentrations, bounded temperatures)
- **Nonlinear constraints**: $c(x_i) \le 0$ (manifold constraints, physical
  coupling)

State constraints are the primary advantage of NMHE over the EKF and UKF,
which cannot enforce hard bounds on state estimates.

In ctrlpp's `nmhe` the box bounds are hard variable bounds on every state in
the window, and the nonlinear constraint is `path_constraint_cb`, `NC` rows
$g(x_i) \le 0$ imposed on every state; with `soft_constraints` it is relaxed
by non-negative slacks under the linear penalty `path_penalty`. A hard bound
$\lvert z_i - h(x_i) \rvert \le r$ on the measurement residuals is also
available. Constraints coupling more than one state are not supported.

## Computational Considerations

### Solution Methods

NMHE is typically solved using [2, Sec. 5, pp. 405--410]:

1. **Sequential Quadratic Programming (SQP)**: iterates over QP subproblems,
   exploiting the structure of the NLP
2. **Interior-point methods**: handle inequality constraints via barrier
   functions

ctrlpp ships SQP-type backends: NLopt's SLSQP (the default of
`nlopt_solver`) and argmin's SQP policies, besides argmin's L-BFGS-B and
augmented Lagrangian. No interior-point backend ships.

### Warm-Starting

Between consecutive time steps, the estimation window shifts by one sample.
The previous solution provides an excellent initial guess for the new NLP:
shift the state trajectory, drop the oldest sample, and append the new
measurement [2, Sec. 5.1, pp. 406--408]. ctrlpp's `nmhe` does this, filling
the newest state with the companion EKF's estimate.

### Real-Time Iteration

As in NMPC, the real-time iteration (RTI) scheme can be applied to NMHE:
perform only a single SQP iteration per time step, relying on the quality
of the warm start. This provides bounded and predictable computation time
[2, Sec. 5.2, pp. 408--410]. ctrlpp does not implement it: every update solves
the window NLP until the backend converges or exhausts its budget.

## Comparison with Other Estimators

| Aspect           | Kalman (EKF/UKF)     | Linear MHE           | NMHE                 |
| ---------------- | -------------------- | -------------------- | -------------------- |
| Model            | Linear/nonlinear     | Linear               | Nonlinear            |
| Constraints      | No                   | Box on states        | General nonlinear    |
| Computation      | $O(n^3)$             | QP                   | NLP                  |
| Optimality       | Optimal (linear)     | Optimal (linear+box) | Local optimum        |
| Noise model      | Gaussian             | Gaussian             | Any (via cost)       |
| Implementation   | Simple               | Moderate             | Complex              |

The table describes the methods. In ctrlpp, `mhe` linearizes a nonlinear
model once per solve, both estimators bound the states and the measurement
residuals, only `nmhe` takes a nonlinear (per-state) constraint, and both use
the fixed quadratic cost above, the Gaussian case.

## References

[1] C. V. Rao, J. B. Rawlings, and D. Q. Mayne, "Constrained State Estimation
for Nonlinear Discrete-Time Systems: Stability and Moving Horizon
Approximations," IEEE Transactions on Automatic Control, vol. 48, no. 2,
pp. 246--258, 2003.

[2] M. Diehl, H. J. Ferreau, and N. Haverbeke, "Efficient Numerical Methods
for Nonlinear MPC and Moving Horizon Estimation," in Nonlinear Model
Predictive Control, LNCIS vol. 384, Springer, 2009, pp. 391--417.

[3] J. B. Rawlings, D. Q. Mayne, and M. Diehl, "Model Predictive Control:
Theory, Computation, and Design," 2nd ed., Nob Hill Publishing, 2017.

[4] D. Simon, "Optimal State Estimation: Kalman, H-Infinity, and Nonlinear
Approaches," Wiley, 2006.
