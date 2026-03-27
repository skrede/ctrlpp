# Nonlinear Moving Horizon Estimation

Nonlinear Moving Horizon Estimation (NMHE) extends the MHE framework to
systems with nonlinear dynamics and measurement models. Instead of solving a
quadratic program, NMHE solves a nonlinear program (NLP) over a sliding
window of recent measurements, enabling state estimation for general nonlinear
constrained systems [1, Sec. 4, pp. 248--258].

NMHE is the estimation dual of nonlinear MPC: both solve NLPs at each time
step, both use receding horizons, and both share the same computational
infrastructure (multiple shooting, SQP solvers).

## NMHE Formulation

For a discrete nonlinear system $x_{i+1} = f(x_i, u_i) + w_i$ with
measurements $z_i = h(x_i) + v_i$, the NMHE problem at time $k$ is
[1, Sec. 4.2, pp. 253--256]:

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
- The measurement residual uses the nonlinear observation $h(x)$ instead of $Cx$
- The resulting optimisation is a nonlinear program, not a QP

## NLP Structure

### Decision Variables

The NLP decision variables are the initial state $x_{k-N}$ and the process
noise sequence $w_{k-N}, \ldots, w_{k-1}$. The remaining states are
determined implicitly by the dynamics constraints
[1, Sec. 4.2, pp. 254--255]:

$$
x_{i+1} = f(x_i, u_i) + w_i, \quad i = k-N, \ldots, k-1
$$

Alternatively, all states can be treated as explicit decision variables
(multiple shooting formulation) with continuity constraints, providing
better convergence properties [2, Sec. 3, pp. 395--400].

### Gradient and Hessian

The NLP gradient requires the Jacobians of $f$ and $h$ with respect to the
state. These can be computed analytically or by finite differences. The Hessian
of the Lagrangian can be approximated using a Gauss-Newton method, which
exploits the least-squares structure of the cost function
[2, Sec. 4, pp. 400--405].

## Arrival Cost for Nonlinear Systems

The arrival cost in NMHE must summarise information from all measurements
before the estimation window. For nonlinear systems, the exact arrival cost
is intractable, so approximations are used [1, Sec. 4.3, pp. 256--258]:

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
prior information but is computationally more expensive.

## State Constraints

NMHE handles state constraints directly in the NLP formulation. Common
constraint types [1, Sec. 4.1, pp. 253--254]:

- **Box constraints**: $x_{\min} \le x_i \le x_{\max}$ (non-negative
  concentrations, bounded temperatures)
- **Nonlinear constraints**: $c(x_i) \le 0$ (manifold constraints, physical
  coupling)

State constraints are the primary advantage of NMHE over the EKF and UKF,
which cannot enforce hard bounds on state estimates.

## Computational Considerations

### Solution Methods

NMHE is typically solved using [2, Sec. 5, pp. 405--410]:

1. **Sequential Quadratic Programming (SQP)**: iterates over QP subproblems,
   exploiting the structure of the NLP
2. **Interior-point methods**: handle inequality constraints via barrier
   functions

### Warm-Starting

Between consecutive time steps, the estimation window shifts by one sample.
The previous solution provides an excellent initial guess for the new NLP:
shift the state trajectory, drop the oldest sample, and append the new
measurement [2, Sec. 5.1, pp. 406--408].

### Real-Time Iteration

As in NMPC, the real-time iteration (RTI) scheme can be applied to NMHE:
perform only a single SQP iteration per time step, relying on the quality
of the warm start. This provides bounded and predictable computation time
[2, Sec. 5.2, pp. 408--410].

## Comparison with Other Estimators

| Aspect           | Kalman (EKF/UKF)     | Linear MHE           | NMHE                 |
| ---------------- | -------------------- | -------------------- | -------------------- |
| Model            | Linear/nonlinear     | Linear               | Nonlinear            |
| Constraints      | No                   | Box on states        | General nonlinear    |
| Computation      | $O(n^3)$             | QP                   | NLP                  |
| Optimality       | Optimal (linear)     | Optimal (linear+box) | Local optimum        |
| Noise model      | Gaussian             | Gaussian             | Any (via cost)       |
| Implementation   | Simple               | Moderate             | Complex              |

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
