# Moving Horizon Estimation

Moving Horizon Estimation (MHE) is the estimation dual of Model Predictive
Control. Instead of optimizing future inputs over a prediction horizon, MHE
optimizes past state estimates over a fixed window of recent measurements. It
provides optimal state estimation for constrained systems where the Kalman
filter cannot enforce physical bounds on states
[1, Sec. 1, p. 246; Sec. 4, pp. 253--254].

Rao, Rawlings, and Mayne [1] gave a general theory of constrained MHE, with
sufficient conditions for its stability. MHE shares the same
receding-horizon philosophy as MPC: solve a finite-window optimization at each
step, use only the most recent estimate, and slide the window forward.

## The MHE Cost Function

At time $k$, MHE estimates the state trajectory over the most recent $N$ time
steps. For a linear system $x_{i+1} = A x_i + B u_i + w_i$ with measurements
$z_i = H x_i + v_i$, the decision variables are the initial state $x_{k-N}$
and the process noise sequence $\{w_{k-N}, \ldots, w_{k-1}\}$
[1, Sec. 1, p. 248]:

$$
\min_{x_{k-N}, \, w_{k-N:k-1}} \;
\lVert x_{k-N} - \bar{x}_{k-N} \rVert^2_{P_{\text{arr}}^{-1}}
+ \sum_{i=k-N}^{k-1} \lVert w_i \rVert^2_{Q^{-1}}
+ \sum_{i=k-N}^{k} \lVert z_i - H x_i \rVert^2_{R^{-1}}
$$

$$
\text{subject to} \quad x_{i+1} = A x_i + B u_i + w_i
$$

$$
x_{\min} \le x_i \le x_{\max}
$$

where the weighted norm is $\lVert a \rVert^2_M = a^\top M \, a$.

ctrlpp's `mhe` solves the same problem with every state $x_{k-N}, \ldots, x_k$
as a decision variable and the process noise eliminated: $w_i$ is replaced by
the residual $x_{i+1} - A x_i - B u_i$ inside the cost, so the dynamics impose
no equality constraint. Its cost carries a factor $\tfrac{1}{2}$ on every
term and a scalar `arrival_cost_weight` $\omega$ on the arrival term,
$\tfrac{\omega}{2} \lVert x_{k-N} - \bar{x}_{k-N} \rVert^2_{P_{\text{arr}}^{-1}}$.
The state bounds can be softened with non-negative slacks under a linear
penalty, and a hard bound $\lvert z_i - H x_i \rvert \le r$ on the
measurement residuals can be added.

The three terms balance:
1. **Arrival cost**: consistency with prior information before the window
2. **Process noise**: penalizing deviations from the model
3. **Measurement fit**: consistency with observations

## Arrival Cost

The first term $\lVert x_{k-N} - \bar{x}_{k-N} \rVert^2_{P_{\text{arr}}^{-1}}$
is the arrival cost. It summarizes all information from measurements before
the estimation window into a single quadratic penalty
[1, Sec. 1, pp. 247--248].

As the window slides forward, the arrival cost must be updated:

- **Prior mean** $\bar{x}_{k-N}$: the best estimate of the state at the
  window boundary, typically from a companion filter
- **Arrival cost weighting** $P_{\text{arr}}$: encodes the uncertainty of
  that estimate

Because $z_{k-N}$ is inside the window, the prior must summarize only the
measurements before it: it is the filter's predicted estimate
$\hat{x}_{k-N \mid k-N-1}$ and covariance $P_{k-N \mid k-N-1}$, taken before
the correction with $z_{k-N}$, so that measurement is not counted twice.

### Exact Update (Linear Systems)

For linear systems without constraints, the exact arrival cost equals the
Kalman filter covariance: $\bar{x}_{k-N} = \hat{x}_{k-N}^{\text{KF}}$ and
$P_{\text{arr}} = P_{k-N}^{\text{KF}}$ [1, Sec. 1, Eqs. (2)--(3), p. 248], both the
predicted quantities at $k-N$.

### Approximate Update (Nonlinear Systems)

For nonlinear systems, an EKF approximation provides the arrival cost. A
companion EKF runs alongside the MHE, and its state estimate and covariance
are used as $\bar{x}$ and $P_{\text{arr}}$ at each window shift
[2, Sec. 4.3, pp. 403--405].

ctrlpp's `mhe` and `nmhe` both embed an `ekf` for this, also for a linear
model, where the EKF's recursion is the Kalman filter's. The EKF runs on its
own: it is corrected with every measurement and never reset to the window
estimate.

## Relationship to the Kalman Filter

For linear systems without state constraints, MHE with infinite window length
recovers the Kalman filter solution exactly. With finite window, MHE
approximates the Kalman filter but can additionally enforce
[1, Sec. 1, pp. 247--248; Sec. 4, p. 253]:

$$
x_{\min} \le x_i \le x_{\max}
$$

This is the key advantage of MHE: state constraints are handled naturally
within the optimization. Physical constraints such as non-negative
concentrations or bounded temperatures can be directly enforced; ctrlpp's
`mhe` takes element-wise bounds on the states and bounds on the measurement
residuals.

With the exact arrival cost above, the finite window loses nothing while no
bound is active: the window's newest estimate is the Kalman filter's estimate
at any window length. ctrlpp's `mhe` uses the Kalman filter's prior as its
arrival cost, and its test suite checks that its newest estimate matches a
Kalman filter's on an unconstrained linear plant.

## QP Formulation

For linear dynamics, the MHE problem is a quadratic program. The structure
is similar to linear MPC with the time direction reversed
[3, Sec. 4.2.2, pp. 281--282]. ctrlpp poses it in the OSQP form its QP backends
take:

$$
\min_z \; \frac{1}{2} z^\top P z + q^\top z
\quad \text{subject to} \quad l \le A_c z \le u
$$

where $z$ collects the states $x_{k-N}, \ldots, x_k$ and, when the state bounds
are softened, their lower and upper slacks. Since the dynamics are in the
cost, there are no equality rows; $A_c$ holds only the state-bound and
residual-bound rows. $P$ is block tridiagonal over the states, coupling each
state only to its neighbors, and the same QP backends as linear MPC solve it.

ctrlpp's `mhe` accepts any dynamics and measurement callables, but it builds
the QP from a single linearization of them per solve, at the window's middle
sample, without the linearization's constant terms: the model it estimates
with is $x_{i+1} = A x_i + B u_i$, $z_i = H x_i$. It is exact for linear models
without an offset; see the [mhe limitations](../api/mpc/mhe.md#limitations).
For a nonlinear model, [NMHE](nmhe.md) keeps the model nonlinear.

## Convergence and Stability

Rao, Rawlings, and Mayne [1, Sec. 3, pp. 250--253] proved that MHE is an
asymptotically stable observer under the following conditions:

1. The system is uniformly observable (their Definition 2.4)
2. The approximate arrival cost meets their conditions C1 and C2: it is
   bounded above by a K-function of the distance from the estimate, and it
   adds no information the data do not hold. The Kalman filter's arrival cost
   meets them for a linear model with convex constraints; for a nonlinear
   model an approximation generally does not, and the paper rescales it online
3. The estimation window is at least the observability horizon, $N \ge N_o$

Without disturbances the estimation error converges to zero; with bounded
disturbances it stays bounded. In practice, window lengths of 10--50 steps are
typical, depending on the system dynamics and measurement rate.

## Tuning Guidelines

- **Window length $N$**: longer windows use more data but increase
  computation. Start with $N = 10$--$20$.
- **$Q$ and $R$**: same interpretation as in the Kalman filter. $Q$
  reflects model uncertainty, $R$ reflects measurement noise. The estimator
  weights with their inverses, so both must be invertible.
- **State bounds**: should reflect true physical constraints, not
  artificially tight bounds (which slow convergence).

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
