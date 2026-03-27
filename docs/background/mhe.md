# MHE Theory

Moving Horizon Estimation (MHE) is the estimation dual of Model Predictive
Control. Instead of optimising future inputs over a prediction horizon, MHE
optimises past state estimates over a fixed window of recent measurements.
It handles state constraints and non-Gaussian noise naturally.

## Moving Horizon Cost Function

MHE estimates the state trajectory over the most recent $N$ time steps by
solving an optimisation problem. At time $k$ the decision variables are the
initial state $x_{k-N}$ and the process noise sequence
$\{w_{k-N}, \ldots, w_{k-1}\}$:

$$
\min_{x_{k-N}, \, w_{k-N:k-1}} \;
\lVert x_{k-N} - \bar{x}_{k-N} \rVert^2_{P_{\text{arr}}^{-1}}
+ \sum_{i=k-N}^{k-1} \lVert w_i \rVert^2_{Q^{-1}}
+ \sum_{i=k-N}^{k} \lVert v_i \rVert^2_{R^{-1}}
$$

$$
\text{subject to} \quad x_{i+1} = f(x_i, u_i) + w_i
$$

$$
z_i = h(x_i) + v_i
$$

$$
x_i \in \mathcal{X} \quad \text{(state constraints)}
$$

where the weighted norm is $\lVert a \rVert^2_M = a^\top M \, a$.

## Arrival Cost

The first term $\lVert x_{k-N} - \bar{x}_{k-N} \rVert^2_{P_{\text{arr}}^{-1}}$
is the *arrival cost*. It summarises all information from measurements before
the estimation window into a single quadratic penalty.

As the window slides forward, the arrival cost must be updated to preserve
historical information. The prior mean $\bar{x}_{k-N}$ and weighting
$P_{\text{arr}}$ are updated at each step:

- **Linear systems**: the exact arrival cost is the Kalman filter covariance,
  so $\bar{x}_{k-N} = \hat{x}_{k-N}^{\text{KF}}$ and
  $P_{\text{arr}} = P_{k-N}^{\text{KF}}$

- **Nonlinear systems**: an EKF approximation is commonly used. ctrlpp
  auto-constructs a companion EKF that runs alongside the MHE to provide
  the arrival cost update.

## Relationship to Kalman Filter

For linear systems without constraints, MHE with infinite window length
recovers the Kalman filter exactly. With finite window, MHE approximates the
Kalman filter but can additionally enforce:

$$
x_{\min} \le x_i \le x_{\max}
$$

This is the key advantage of MHE over the Kalman filter -- state constraints
are handled naturally within the optimisation.

## QP/NLP Formulation

- **Linear MHE**: linear dynamics $x_{i+1} = A x_i + B u_i + w_i$ yield a
  quadratic program. ctrlpp uses OSQP as the default QP solver.

- **Nonlinear MHE**: nonlinear dynamics $x_{i+1} = f(x_i, u_i) + w_i$ yield
  a nonlinear program. ctrlpp uses NLopt as the default NLP solver.

Both formulations support box constraints on states.

## References

- Rao, C. V., Rawlings, J. B., and Mayne, D. Q. (2003). "Constrained State
  Estimation for Nonlinear Discrete-Time Systems: Stability and Moving
  Horizon Approximations." *IEEE Trans. Automatic Control*,
  48(2):246--258. [`rao2003`]
  Establishes the theoretical foundation for MHE: stability, convergence, and
  the role of the arrival cost.

- Diehl, M., Ferreau, H. J., and Haverbeke, N. (2009). "Efficient Numerical
  Methods for Nonlinear MPC and Moving Horizon Estimation." In *Nonlinear
  Model Predictive Control*, LNCIS vol. 384, pp. 391--417. Springer.
  [`diehl2009`]
  Covers computational aspects of MHE including real-time iteration schemes.

- Rawlings, J. B., Mayne, D. Q., and Diehl, M. (2017). *Model Predictive
  Control: Theory, Computation, and Design.* Nob Hill Publishing, 2nd ed.
  [`rawlings2017`]
  Unified treatment of MPC and MHE as dual optimisation problems.

## Related API Pages

- [mhe](../api/mpc/mhe.md) -- linear moving horizon estimation (OSQP)
- [nmhe](../api/mpc/nmhe.md) -- nonlinear moving horizon estimation (NLopt)
- [MPC Theory](mpc.md) -- the control dual of MHE
