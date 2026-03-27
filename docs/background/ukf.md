# Unscented Kalman Filter

The Unscented Kalman Filter (UKF) is a derivative-free nonlinear state
estimator that propagates a set of carefully chosen sample points (sigma
points) through the true nonlinear functions, then reconstructs the
transformed mean and covariance from the propagated points. This avoids
the Jacobian computation required by the EKF and captures second-order
statistics accurately for any nonlinearity [1, Ch. 14, pp. 433--459].

The UKF was introduced by Julier and Uhlmann [2] and later refined by
Wan and Van der Merwe [3] with the scaled unscented transform. It has
become a standard alternative to the EKF for moderately nonlinear systems
where Jacobians are expensive or unavailable.

## The Unscented Transform

The unscented transform (UT) is the core mechanism of the UKF. Given a
random variable $x \in \mathbb{R}^n$ with mean $\hat{x}$ and covariance $P$,
the UT approximates the statistics of $y = g(x)$ for an arbitrary nonlinear
function $g$ [2, Sec. III, pp. 404--408].

### Sigma-Point Generation

For an $n$-dimensional state, the scaled UT generates $2n + 1$ sigma points
[3, Sec. 7.2, pp. 226--230]:

$$
\mathcal{X}_0 = \hat{x}
$$

$$
\mathcal{X}_i = \hat{x} + \left(\sqrt{(n + \lambda) \, P}\right)_i, \quad i = 1, \ldots, n
$$

$$
\mathcal{X}_{i+n} = \hat{x} - \left(\sqrt{(n + \lambda) \, P}\right)_i, \quad i = 1, \ldots, n
$$

where $\lambda = \alpha^2 (n + \kappa) - n$ is the composite scaling parameter
and $\left(\sqrt{M}\right)_i$ denotes the $i$-th column of a matrix square
root of $M$ (typically computed via Cholesky or LDLT decomposition).

### Sigma-Point Weights

The mean and covariance weights are [3, Sec. 7.2, pp. 228--229]:

$$
W_0^{(m)} = \frac{\lambda}{n + \lambda}, \qquad
W_0^{(c)} = \frac{\lambda}{n + \lambda} + (1 - \alpha^2 + \beta)
$$

$$
W_i^{(m)} = W_i^{(c)} = \frac{1}{2(n + \lambda)}, \quad i = 1, \ldots, 2n
$$

### Scaling Parameters

The three parameters control the sigma-point distribution
[3, Sec. 7.2, pp. 229--230]:

- **$\alpha$**: determines the spread of sigma points around the mean.
  Small values ($10^{-4}$ to $1$) keep points close to the mean.
- **$\beta$**: incorporates prior knowledge of the distribution.
  $\beta = 2$ is optimal for Gaussian distributions.
- **$\kappa$**: secondary scaling parameter. A common choice is
  $\kappa = 3 - n$ to match the fourth moment of the Gaussian
  [2, Sec. IV, p. 410].

### Transform Procedure

Given sigma points $\mathcal{X}_i$ and a nonlinear function $g$
[3, Sec. 7.3, pp. 231--234]:

1. **Propagate**: $\mathcal{Y}_i = g(\mathcal{X}_i)$ for $i = 0, \ldots, 2n$
2. **Recover mean**: $\bar{y} = \sum_{i=0}^{2n} W_i^{(m)} \mathcal{Y}_i$
3. **Recover covariance**: $P_{yy} = \sum_{i=0}^{2n} W_i^{(c)} (\mathcal{Y}_i - \bar{y})(\mathcal{Y}_i - \bar{y})^\top$
4. **Cross-covariance**: $P_{xy} = \sum_{i=0}^{2n} W_i^{(c)} (\mathcal{X}_i - \hat{x})(\mathcal{Y}_i - \bar{y})^\top$

This captures the mean and covariance accurately to second order for any
nonlinearity, and to third order for Gaussian inputs [2, Sec. V, p. 412].

## UKF Algorithm

The UKF applies the unscented transform twice per time step: once for
prediction and once for the measurement update
[1, Sec. 14.3, pp. 442--448].

### Prediction

1. Generate sigma points from $(\hat{x}_{k-1|k-1}, P_{k-1|k-1})$
2. Propagate through the process model:
   $\mathcal{X}_{k|k-1}^{(i)} = f(\mathcal{X}_{k-1|k-1}^{(i)}, u_{k-1})$
3. Compute the predicted mean and covariance:

$$
\hat{x}_{k|k-1} = \sum_i W_i^{(m)} \mathcal{X}_{k|k-1}^{(i)}
$$

$$
P_{k|k-1} = \sum_i W_i^{(c)} (\mathcal{X}_{k|k-1}^{(i)} - \hat{x}_{k|k-1})(\mathcal{X}_{k|k-1}^{(i)} - \hat{x}_{k|k-1})^\top + Q
$$

### Update

1. Regenerate sigma points from $(\hat{x}_{k|k-1}, P_{k|k-1})$
2. Transform through the measurement model:
   $\mathcal{Z}_{k}^{(i)} = h(\mathcal{X}_{k|k-1}^{(i)})$
3. Compute predicted measurement and covariances:

$$
\hat{z}_k = \sum_i W_i^{(m)} \mathcal{Z}_k^{(i)}
$$

$$
P_{zz} = \sum_i W_i^{(c)} (\mathcal{Z}_k^{(i)} - \hat{z}_k)(\mathcal{Z}_k^{(i)} - \hat{z}_k)^\top + R
$$

$$
P_{xz} = \sum_i W_i^{(c)} (\mathcal{X}_{k|k-1}^{(i)} - \hat{x}_{k|k-1})(\mathcal{Z}_k^{(i)} - \hat{z}_k)^\top
$$

4. Compute the Kalman gain and update:

$$
K_k = P_{xz} \, P_{zz}^{-1}
$$

$$
\hat{x}_{k|k} = \hat{x}_{k|k-1} + K_k (z_k - \hat{z}_k)
$$

$$
P_{k|k} = P_{k|k-1} - K_k \, P_{zz} \, K_k^\top
$$

## Sigma-Point Strategies

Two main strategies exist for generating sigma points:

### Van der Merwe (Scaled)

Uses three parameters $(\alpha, \beta, \kappa)$ for maximum flexibility. The
scaling allows placing sigma points arbitrarily close to or far from the mean.
This is the most widely used variant [3, Sec. 7.2].

### Julier (Original)

The original unscented transform uses a single parameter $\kappa$, with
$\alpha = 1$ and $\beta = 0$ implicitly [2, Sec. III]. Simpler but less
flexible. The zeroth sigma point has weight $\kappa / (n + \kappa)$, which
can be negative for large state dimensions if $\kappa < 0$.

## EKF vs UKF Comparison

| Criterion              | EKF                        | UKF                        |
| ---------------------- | -------------------------- | -------------------------- |
| Nonlinearity           | Mild                       | Moderate to strong         |
| Jacobians required     | Yes                        | No                         |
| Accuracy               | First order                | Second order               |
| Computational cost     | $O(n^2)$ per Jacobian      | $O(n^3)$ for Cholesky      |
| State dimension        | Any                        | Moderate ($2n+1$ points)   |
| Implementation effort  | Jacobians needed           | Only dynamics function     |

For mildly nonlinear systems where analytical Jacobians are available, the EKF
is typically sufficient and more efficient. For strongly nonlinear systems or
when Jacobians are expensive to derive, the UKF is preferred
[1, Sec. 14.5, pp. 455--459].

## References

[1] D. Simon, "Optimal State Estimation: Kalman, H-Infinity, and Nonlinear
Approaches," Wiley, 2006.

[2] S. J. Julier and J. K. Uhlmann, "Unscented Filtering and Nonlinear
Estimation," Proceedings of the IEEE, vol. 92, no. 3, pp. 401--422, 2004.

[3] E. A. Wan and R. van der Merwe, "The Unscented Kalman Filter," in Kalman
Filtering and Neural Networks, S. Haykin, Ed., Wiley, 2001, Ch. 7,
pp. 221--280.

[4] R. van der Merwe, "Sigma-Point Kalman Filters for Probabilistic Inference
in Dynamic State-Space Models," Ph.D. dissertation, Oregon Health & Science
University, 2004.
