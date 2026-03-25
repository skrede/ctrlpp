# EKF and UKF Theory

When the system dynamics or measurement model are nonlinear, the standard
Kalman filter no longer applies. The Extended Kalman Filter (EKF) and
Unscented Kalman Filter (UKF) extend the Kalman framework to nonlinear
systems using different approximation strategies.

## Nonlinear System Model

Both filters address the general discrete nonlinear system:

$$
x_{k+1} = f(x_k, u_k) + w_k, \qquad w_k \sim \mathcal{N}(0, Q)
$$

$$
z_k = h(x_k) + v_k, \qquad v_k \sim \mathcal{N}(0, R)
$$

where $f$ is the process model and $h$ is the measurement model.

## EKF: Linearisation Approach

The EKF applies the standard Kalman equations to a first-order Taylor
expansion of the nonlinear functions. At each step the Jacobians are evaluated
at the current estimate:

$$
F_k = \left. \frac{\partial f}{\partial x} \right|_{\hat{x}_{k|k}}, \qquad
H_k = \left. \frac{\partial h}{\partial x} \right|_{\hat{x}_{k|k-1}}
$$

### EKF Prediction

$$
\hat{x}_{k+1|k} = f(\hat{x}_{k|k}, u_k)
$$

$$
P_{k+1|k} = F_k \, P_{k|k} \, F_k^\top + Q
$$

### EKF Update

$$
K_k = P_{k|k-1} \, H_k^\top \, (H_k \, P_{k|k-1} \, H_k^\top + R)^{-1}
$$

$$
\hat{x}_{k|k} = \hat{x}_{k|k-1} + K_k \bigl(z_k - h(\hat{x}_{k|k-1})\bigr)
$$

$$
P_{k|k} = (I - K_k \, H_k) \, P_{k|k-1}
$$

**Strengths**: simple, well-understood, efficient for mildly nonlinear systems.

**Limitations**: first-order approximation can diverge for strongly nonlinear
systems. Requires Jacobian computation (analytical or numerical). ctrlpp
supports both via `if constexpr` dispatch -- if the user provides Jacobian
functions they are used; otherwise finite-difference approximations are
computed automatically.

## UKF: Sigma-Point Transform

The UKF avoids linearisation entirely. Instead, it propagates a set of
carefully chosen sample points (sigma points) through the true nonlinear
functions and reconstructs the mean and covariance from the transformed
points.

### Sigma-Point Generation

For an $n$-dimensional state with mean $\hat{x}$ and covariance $P$, the
scaled unscented transform generates $2n + 1$ sigma points:

$$
\mathcal{X}_0 = \hat{x}
$$

$$
\mathcal{X}_i = \hat{x} + \left(\sqrt{(n + \lambda) \, P}\right)_i, \quad i = 1, \ldots, n
$$

$$
\mathcal{X}_{i+n} = \hat{x} - \left(\sqrt{(n + \lambda) \, P}\right)_i, \quad i = 1, \ldots, n
$$

where $\lambda = \alpha^2 (n + \kappa) - n$ is the scaling parameter and
$\left(\sqrt{M}\right)_i$ denotes the $i$-th column of the matrix square root
(computed via LDLT decomposition in ctrlpp).

### Sigma-Point Weights

The mean and covariance weights are:

$$
W_0^{(m)} = \frac{\lambda}{n + \lambda}, \qquad
W_0^{(c)} = \frac{\lambda}{n + \lambda} + (1 - \alpha^2 + \beta)
$$

$$
W_i^{(m)} = W_i^{(c)} = \frac{1}{2(n + \lambda)}, \quad i = 1, \ldots, 2n
$$

The parameter $\beta = 2$ is optimal for Gaussian distributions.

### Unscented Transform

Given a nonlinear function $g$ and sigma points $\mathcal{X}_i$ with
weights $W_i$:

1. Propagate each sigma point: $\mathcal{Y}_i = g(\mathcal{X}_i)$
2. Recover the transformed mean: $\bar{y} = \sum_i W_i^{(m)} \mathcal{Y}_i$
3. Recover the transformed covariance: $P_{yy} = \sum_i W_i^{(c)} (\mathcal{Y}_i - \bar{y})(\mathcal{Y}_i - \bar{y})^\top$

This captures mean and covariance accurately to second order for any
nonlinearity, without computing Jacobians.

### Sigma-Point Strategies

ctrlpp implements two sigma-point strategies as swappable template parameters:

- **Van der Merwe (scaled)**: uses parameters $\alpha$, $\beta$, $\kappa$ to
  control the spread of sigma points. This is the default strategy.

- **Julier**: the original unscented transform. Uses a single $\kappa$
  parameter. Simpler but less flexible than the scaled variant.

## When to Prefer EKF vs UKF

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

- Julier, S. J. and Uhlmann, J. K. (2004). "Unscented Filtering and
  Nonlinear Estimation." *Proceedings of the IEEE*, 92(3):401--422.
  [`julier2004`]
  The definitive paper on the unscented transform and its application to
  state estimation.

- Wan, E. A. and van der Merwe, R. (2001). "The Unscented Kalman Filter."
  In *Kalman Filtering and Neural Networks*, Ch. 7, pp. 221--280. Wiley.
  [`wan2001`]
  Introduces the scaled unscented transform with $\alpha$/$\beta$/$\kappa$
  parameters.

- Van der Merwe, R. (2004). *Sigma-Point Kalman Filters for Probabilistic
  Inference in Dynamic State-Space Models.* PhD thesis, Oregon Health &
  Science University. [`vandermerwe2004`]
  Comprehensive treatment of sigma-point methods including square-root
  variants.

## Related API Pages

- [ekf](../estimation/ekf.md) -- Extended Kalman filter with analytical/
  numerical Jacobian dispatch
- [ukf](../estimation/ukf.md) -- Unscented Kalman filter with swappable
  sigma-point strategies
- [Kalman Theory](kalman-theory.md) -- linear estimation foundation
