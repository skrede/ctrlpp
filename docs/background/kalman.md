# Kalman Filter Theory

The Kalman filter is the optimal linear estimator for systems with Gaussian
noise. It recursively estimates the state of a dynamic system from a sequence
of noisy measurements, producing the minimum-variance unbiased estimate at
each time step.

## State Estimation Problem

Given a discrete linear time-invariant system:

$$
x_{k+1} = A \, x_k + B \, u_k + w_k
$$

$$
z_k = C \, x_k + v_k
$$

where $w_k \sim \mathcal{N}(0, Q)$ is process noise and
$v_k \sim \mathcal{N}(0, R)$ is measurement noise, the Kalman filter computes
the state estimate $\hat{x}_k$ that minimises the expected squared estimation
error $\mathbb{E}\bigl[\lVert x_k - \hat{x}_k \rVert^2\bigr]$.

## Prediction Step

The prediction propagates the state estimate and error covariance forward:

$$
\hat{x}_{k|k-1} = A \, \hat{x}_{k-1|k-1} + B \, u_{k-1}
$$

$$
P_{k|k-1} = A \, P_{k-1|k-1} \, A^\top + Q
$$

The covariance $P_{k|k-1}$ grows during prediction, reflecting increased
uncertainty about the state.

## Update (Correction) Step

When a measurement $z_k$ arrives, the filter computes the innovation
(measurement residual):

$$
\tilde{y}_k = z_k - C \, \hat{x}_{k|k-1}
$$

The innovation covariance is:

$$
S_k = C \, P_{k|k-1} \, C^\top + R
$$

The Kalman gain optimally blends prediction and measurement:

$$
K_k = P_{k|k-1} \, C^\top \, S_k^{-1}
$$

The corrected state estimate and covariance are:

$$
\hat{x}_{k|k} = \hat{x}_{k|k-1} + K_k \, \tilde{y}_k
$$

$$
P_{k|k} = (I - K_k \, C) \, P_{k|k-1}
$$

The Joseph form of the covariance update is numerically more stable:

$$
P_{k|k} = (I - K_k C) \, P_{k|k-1} \, (I - K_k C)^\top + K_k \, R \, K_k^\top
$$

## Optimality

For linear systems with Gaussian noise, the Kalman filter is:

- **Minimum variance**: no other linear unbiased estimator achieves smaller
  estimation error covariance
- **Maximum likelihood**: the estimate is the most probable state given all
  measurements
- **Recursive**: processes measurements one at a time without storing history

## Covariance Propagation

The error covariance $P$ encodes the estimator's confidence. It depends only
on the system matrices $(A, C)$ and noise covariances $(Q, R)$, not on the
actual measurements. This means the Kalman gain sequence can be precomputed
for time-invariant systems, converging to the steady-state gain
$K_\infty = P_\infty C^\top (C P_\infty C^\top + R)^{-1}$ where $P_\infty$
satisfies the discrete algebraic Riccati equation.

## References

- Kalman, R. E. (1960). "A New Approach to Linear Filtering and Prediction
  Problems." *Journal of Basic Engineering*, 82(1):35--45. [`kalman1960`]
  The foundational paper establishing recursive optimal estimation for linear
  systems.

- Simon, D. (2006). *Optimal State Estimation: Kalman, H-Infinity, and
  Nonlinear Approaches.* Wiley. [`simon2006`]
  Comprehensive reference covering the Kalman filter, its extensions, and
  connections to other estimation frameworks.

- Kailath, T., Sayed, A. H., and Hassibi, B. (2000). *Linear Estimation.*
  Prentice Hall. [`kailath2000`]
  Rigorous treatment of linear estimation theory including square-root and
  information filter variants.

## Related API Pages

- [kalman_filter](../api/estimation/kalman.md) -- linear Kalman filter
  implementation
- [luenberger_observer](../api/estimation/luenberger.md) -- deterministic
  observer (Kalman filter without stochastic tuning)
- [Your First Estimator](../guides/intro/your-first-estimator.md) --
  introductory tutorial
- [Observer-Controller Guide](../guides/estimation/observer-controller.md) --
  composing observers with controllers
