# EKF and UKF Theory

When the system dynamics or measurement model are nonlinear, the standard
Kalman filter no longer applies. The Extended Kalman Filter (EKF) and
Unscented Kalman Filter (UKF) extend the Kalman framework to nonlinear
systems using different approximation strategies.

## Key Concepts

### EKF: Linearisation Approach

The EKF applies the standard Kalman equations to a first-order Taylor
expansion of the nonlinear functions:

```
x[k+1] = f(x[k], u[k]) + w[k]
z[k]   = h(x[k]) + v[k]
```

At each step, the Jacobians F = df/dx and H = dh/dx are evaluated at the
current estimate. The prediction and update use these linearised matrices in
place of A and C.

**Strengths**: simple, well-understood, efficient for mildly nonlinear systems.

**Limitations**: first-order approximation can diverge for strongly nonlinear
systems. Requires Jacobian computation (analytical or numerical).

### UKF: Sigma-Point Transform

The UKF avoids linearisation entirely. Instead, it propagates a set of
carefully chosen sample points (sigma points) through the true nonlinear
functions and reconstructs the mean and covariance from the transformed
points.

For an n-dimensional state, the UKF generates 2n+1 sigma points that capture
the mean and covariance of the current estimate. Each point is propagated
through f() or h(), and the statistics are recovered via weighted averaging.

**Strengths**: captures second-order effects, no Jacobians needed, better
accuracy for strongly nonlinear systems.

**Limitations**: higher computational cost (2n+1 function evaluations per
step), requires careful sigma-point parameter tuning.

### Sigma-Point Strategies

ctrlpp implements two sigma-point strategies as swappable template parameters:

- **Van der Merwe (scaled)**: uses parameters alpha, beta, kappa to control
  the spread of sigma points. Alpha controls the spread around the mean;
  beta incorporates prior distribution knowledge (beta=2 is optimal for
  Gaussian). This is the default strategy.

- **Julier**: the original unscented transform. Uses a single kappa parameter.
  Simpler but less flexible than the scaled variant.

### When to Prefer EKF vs UKF

| Criterion              | EKF                        | UKF                        |
| ---------------------- | -------------------------- | -------------------------- |
| Nonlinearity           | Mild                       | Moderate to strong         |
| Jacobians available    | Yes (analytical preferred) | Not needed                 |
| Computational budget   | Tighter                    | More generous              |
| State dimension        | Any                        | Moderate (sigma pt. count) |
| Implementation effort  | Needs Jacobians            | Just dynamics function     |

For mildly nonlinear systems where analytical Jacobians are available, the EKF
is typically sufficient and more efficient. For strongly nonlinear systems or
when Jacobians are expensive to compute, the UKF is the better choice.

## References

- **Julier, S. J. and Uhlmann, J. K.** "Unscented Filtering and Nonlinear
  Estimation." *Proceedings of the IEEE*, 92(3):401--422, 2004.
  DOI: 10.1109/JPROC.2003.823141.
  The definitive paper on the unscented transform and its application to
  state estimation.

- **Wan, E. A. and van der Merwe, R.** "The Unscented Kalman Filter." In
  *Kalman Filtering and Neural Networks*, Ch. 7, pp. 221--280. Wiley, 2001.
  DOI: 10.1002/0471221546.ch7.
  Introduces the scaled unscented transform with alpha/beta/kappa parameters.

- **Van der Merwe, R.** *Sigma-Point Kalman Filters for Probabilistic
  Inference in Dynamic State-Space Models.* PhD thesis, Oregon Health &
  Science University, 2004.
  Comprehensive treatment of sigma-point methods including square-root
  variants.

## Related API Pages

- [ekf](../estimation/ekf.md) -- Extended Kalman filter with analytical/
  numerical Jacobian dispatch
- [ukf](../estimation/ukf.md) -- Unscented Kalman filter with swappable
  sigma-point strategies
- [Kalman Theory](kalman-theory.md) -- linear estimation foundation
- [Observer-Controller Guide](../guides/estimation/observer-controller.md) -- composing observers with controllers
