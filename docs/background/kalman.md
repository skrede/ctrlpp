# Kalman Filter

The Kalman filter is the optimal linear estimator for systems with Gaussian
noise. It recursively estimates the state of a dynamic system from a sequence
of noisy measurements, producing the minimum-variance unbiased estimate at
each time step [1, Ch. 5, pp. 123--148].

The filter was introduced by Kalman in 1960 [2] and has since become the
foundation of state estimation in navigation, aerospace, robotics, and
control. Its recursive structure makes it suitable for real-time
implementation, processing each measurement as it arrives without storing
the full measurement history.

## State Estimation Problem

Consider a discrete linear time-invariant (LTI) system
[1, Sec. 5.1, pp. 124--125]:

$$
x_{k+1} = A \, x_k + B \, u_k + w_k
$$

$$
z_k = C \, x_k + v_k
$$

where $x_k \in \mathbb{R}^n$ is the state, $u_k$ is the known control input,
$z_k$ is the measurement, $w_k \sim \mathcal{N}(0, Q)$ is process noise, and
$v_k \sim \mathcal{N}(0, R)$ is measurement noise. The noise sequences are
assumed white, zero-mean, and mutually uncorrelated.

The goal is to compute the estimate $\hat{x}_{k|k}$ that minimises the
expected squared estimation error
$\mathbb{E}\bigl[\lVert x_k - \hat{x}_k \rVert^2\bigr]$.

## Prediction Step

The prediction propagates the state estimate and error covariance forward
through the system dynamics [1, Sec. 5.2, pp. 128--130]:

$$
\hat{x}_{k|k-1} = A \, \hat{x}_{k-1|k-1} + B \, u_{k-1}
$$

$$
P_{k|k-1} = A \, P_{k-1|k-1} \, A^\top + Q
$$

The predicted covariance $P_{k|k-1}$ grows during prediction, reflecting
increased uncertainty about the state between measurements.

## Update (Correction) Step

When a measurement $z_k$ arrives, the filter computes the innovation
(measurement residual) [1, Sec. 5.2, pp. 130--133]:

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

## Joseph Form

The standard covariance update $P_{k|k} = (I - K_k C) P_{k|k-1}$ is
numerically sensitive. The Joseph form provides guaranteed symmetry and
positive semi-definiteness [1, Sec. 5.4, pp. 139--140]:

$$
P_{k|k} = (I - K_k C) \, P_{k|k-1} \, (I - K_k C)^\top + K_k \, R \, K_k^\top
$$

This form is more expensive but essential for long-running filters where
numerical drift can accumulate.

## Optimality Properties

For linear systems with Gaussian noise, the Kalman filter is
[1, Sec. 5.3, pp. 134--138]:

- **Minimum variance**: no other linear unbiased estimator achieves smaller
  estimation error covariance.
- **Maximum likelihood**: the estimate is the most probable state given all
  measurements.
- **Recursive**: processes measurements one at a time without storing history.
- **Unbiased**: the expected value of the estimation error is zero.

For non-Gaussian noise, the Kalman filter remains the best linear unbiased
estimator (BLUE), though nonlinear estimators may perform better.

## Discrete Algebraic Riccati Equation

For time-invariant systems, the error covariance $P$ converges to a
steady-state value $P_\infty$ that satisfies the discrete algebraic Riccati
equation (DARE) [1, Sec. 5.5, pp. 141--145]:

$$
P_\infty = A \, P_\infty \, A^\top - A \, P_\infty \, C^\top (C \, P_\infty \, C^\top + R)^{-1} C \, P_\infty \, A^\top + Q
$$

The steady-state Kalman gain is then:

$$
K_\infty = P_\infty \, C^\top (C \, P_\infty \, C^\top + R)^{-1}
$$

This allows precomputing the gain, reducing the online computation to a
matrix-vector multiply. Convergence requires the system to be detectable
(all unstable modes are observable) [1, Sec. 5.5, p. 143].

## Noise Tuning

The performance of the Kalman filter depends critically on the choice of
process noise covariance $Q$ and measurement noise covariance $R$
[1, Sec. 5.6, pp. 146--148]:

- **$Q$ too small**: the filter trusts the model too much and responds slowly
  to actual disturbances.
- **$Q$ too large**: the filter trusts measurements too much, becoming noisy.
- **$R$ too small**: the filter over-reacts to measurement noise.
- **$R$ too large**: the filter ignores measurements and relies on prediction.

The ratio $Q/R$ determines the bandwidth of the estimator. In practice, $R$
can often be measured directly from sensor data sheets, while $Q$ requires
tuning based on the expected magnitude of model uncertainty and disturbances.

## Observability

The Kalman filter converges only if the system is observable -- that is, the
state can be uniquely determined from a finite sequence of measurements. The
observability matrix [3, Sec. 12.3, pp. 683--690]:

$$
\mathcal{O} = \begin{bmatrix} C \\ CA \\ CA^2 \\ \vdots \\ CA^{n-1} \end{bmatrix}
$$

must have rank $n$ (full column rank). If the system is not fully observable,
the Kalman filter can only estimate the observable subspace.

## References

[1] D. Simon, "Optimal State Estimation: Kalman, H-Infinity, and Nonlinear
Approaches," Wiley, 2006.

[2] R. E. Kalman, "A New Approach to Linear Filtering and Prediction
Problems," Journal of Basic Engineering, vol. 82, no. 1, pp. 35--45, 1960.

[3] N. S. Nise, "Control Systems Engineering," 7th ed., Wiley, 2015.

[4] S. Sarkka and L. Svensson, "Bayesian Filtering and Smoothing," 2nd ed.,
Cambridge University Press, 2023.
