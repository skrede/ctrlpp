# Kalman Filter Theory

The Kalman filter is the optimal linear estimator for systems with Gaussian
noise. It recursively estimates the state of a dynamic system from a sequence
of noisy measurements, producing the minimum-variance unbiased estimate at
each time step.

## Key Concepts

### State Estimation Problem

Given a discrete linear system:

```
x[k+1] = A * x[k] + B * u[k] + w[k]     (process model)
z[k]   = C * x[k] + v[k]                  (measurement model)
```

where w[k] ~ N(0, Q) and v[k] ~ N(0, R) are zero-mean Gaussian noise, the
Kalman filter finds the state estimate x_hat[k] that minimises the expected
squared estimation error.

### Prediction Step

The prediction propagates the state estimate and covariance forward in time:

```
x_hat[k|k-1] = A * x_hat[k-1|k-1] + B * u[k-1]
P[k|k-1]     = A * P[k-1|k-1] * A^T + Q
```

The covariance P grows during prediction, reflecting increased uncertainty.

### Update (Correction) Step

When a measurement arrives, the Kalman gain determines the optimal blend
between prediction and measurement:

```
K[k]         = P[k|k-1] * C^T * (C * P[k|k-1] * C^T + R)^{-1}
x_hat[k|k]   = x_hat[k|k-1] + K[k] * (z[k] - C * x_hat[k|k-1])
P[k|k]       = (I - K[k] * C) * P[k|k-1]
```

The innovation z[k] - C * x_hat[k|k-1] measures how much the measurement
deviates from the prediction. The Kalman gain weights this innovation
according to the relative uncertainties.

### Optimality

For linear systems with Gaussian noise, the Kalman filter is:

- **Minimum variance**: no other linear unbiased estimator achieves smaller
  estimation error covariance
- **Maximum likelihood**: the estimate is the most probable state given all
  measurements
- **Recursive**: processes measurements one at a time without storing
  history

### Covariance Propagation

The error covariance P encodes the estimator's confidence. It depends only on
the system matrices and noise covariances, not on the actual measurements.
This means the Kalman gain sequence can be precomputed for time-invariant
systems.

## References

- **Kalman, R. E.** "A New Approach to Linear Filtering and Prediction
  Problems." *Journal of Basic Engineering*, 82(1):35--45, 1960.
  DOI: 10.1115/1.3662552.
  The foundational paper establishing recursive optimal estimation for linear
  systems.

- **Simon, D.** *Optimal State Estimation: Kalman, H-Infinity, and Nonlinear
  Approaches.* Wiley, 2006. ISBN 978-0-471-70858-2.
  Comprehensive reference covering the Kalman filter, its extensions, and
  connections to other estimation frameworks.

- **Kailath, T., Sayed, A. H., and Hassibi, B.** *Linear Estimation.*
  Prentice Hall, 2000. ISBN 978-0-13-022464-4.
  Rigorous treatment of linear estimation theory including square-root and
  information filter variants.

## Related API Pages

- [kalman_filter](../estimation/kalman.md) -- linear Kalman filter
  implementation
- [luenberger_observer](../estimation/luenberger.md) -- deterministic
  observer (Kalman filter without stochastic tuning)
- [Your First Estimator](../guides/intro/your-first-estimator.md) --
  introductory tutorial
