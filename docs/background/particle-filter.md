# Particle Filter Theory

Particle filters (sequential Monte Carlo methods) represent the posterior
distribution of the state using a set of weighted samples (particles). Unlike
the Kalman family, particle filters make no assumptions about linearity or
Gaussianity, making them applicable to arbitrary distributions and strongly
nonlinear systems.

## Importance Sampling

The posterior distribution $p(x_k \mid z_{1:k})$ is approximated by $N_p$
weighted particles:

$$
p(x_k \mid z_{1:k}) \approx \sum_{i=1}^{N_p} w_k^{(i)} \, \delta(x_k - x_k^{(i)})
$$

where $x_k^{(i)}$ is the $i$-th particle, $w_k^{(i)}$ is its normalised
weight, and $\delta$ is the Dirac delta. Each particle represents a
hypothesis about the true state; particles with higher weights are more
consistent with the observations.

## Bootstrap Filter (SIR)

The Sequential Importance Resampling (SIR) or bootstrap filter is the
simplest and most widely used particle filter variant. ctrlpp implements this
algorithm:

**1. Propagate** -- draw each particle through the dynamics model with process
noise:

$$
x_k^{(i)} \sim p(x_k \mid x_{k-1}^{(i)}, u_{k-1})
$$

In practice this means evaluating $x_k^{(i)} = f(x_{k-1}^{(i)}, u_{k-1}) + w_k^{(i)}$
where $w_k^{(i)} \sim \mathcal{N}(0, Q)$.

**2. Weight** -- compute the likelihood of the measurement given each
particle:

$$
\tilde{w}_k^{(i)} = p(z_k \mid x_k^{(i)})
$$

For Gaussian measurement noise $v_k \sim \mathcal{N}(0, R)$:

$$
\tilde{w}_k^{(i)} \propto \exp\!\left(-\frac{1}{2} (z_k - h(x_k^{(i)}))^\top R^{-1} (z_k - h(x_k^{(i)}))\right)
$$

**3. Normalise** -- scale weights to sum to one:

$$
w_k^{(i)} = \frac{\tilde{w}_k^{(i)}}{\sum_{j=1}^{N_p} \tilde{w}_k^{(j)}}
$$

**4. Resample** -- draw $N_p$ particles with replacement, proportional to
weights (when triggered by ESS criterion).

## Effective Sample Size and Resampling

After several steps, most particles accumulate negligible weight while a few
carry all the mass. This *weight degeneracy* is measured by the effective
sample size:

$$
N_{\text{eff}} = \frac{1}{\sum_{i=1}^{N_p} \left(w_k^{(i)}\right)^2}
$$

When $N_{\text{eff}}$ drops below a threshold (typically $N_p / 2$),
resampling is triggered. ctrlpp uses systematic resampling: generate a single
random offset $u_0 \sim \mathcal{U}(0, 1/N_p)$ and take samples at
$u_0 + j/N_p$ for $j = 0, \ldots, N_p - 1$ from the cumulative weight
distribution. This has lower variance than multinomial resampling.

## Roughening

After resampling, many particles may be identical copies. *Roughening* adds
small random perturbations to maintain diversity:

$$
x_k^{(i)} \leftarrow x_k^{(i)} + \epsilon^{(i)}, \qquad \epsilon^{(i)} \sim \mathcal{N}(0, \sigma_r^2 I)
$$

where $\sigma_r$ is proportional to the spread of the particle cloud. This
prevents sample impoverishment and helps the filter explore the state space.

## Convergence

As $N_p \to \infty$, the particle approximation converges to the true
posterior. In practice, the required $N_p$ grows exponentially with state
dimension -- the "curse of dimensionality". Particle filters are most
effective for low-dimensional state spaces (typically fewer than 6--10
dimensions).

## References

- Gordon, N. J., Salmond, D. J., and Smith, A. F. M. (1993). "Novel
  Approach to Nonlinear/Non-Gaussian Bayesian State Estimation." *IEE
  Proceedings F*, 140(2):107--113. [`gordon1993`]
  The foundational paper introducing the bootstrap particle filter.

- Arulampalam, M. S., Maskell, S., Gordon, N., and Clapp, T. (2002). "A
  Tutorial on Particle Filters for Online Nonlinear/Non-Gaussian Bayesian
  Tracking." *IEEE Transactions on Signal Processing*, 50(2):174--188.
  [`arulampalam2002`]
  Comprehensive tutorial covering SIR, auxiliary particle filters, and
  practical implementation considerations.

## Related API Pages

- [particle_filter](../api/estimation/particle-filter.md) -- bootstrap SIR
  particle filter with compile-time particle count and ESS-adaptive
  resampling
- [EKF/UKF Theory](ekf-ukf.md) -- alternative nonlinear estimators
- [Observer-Controller Guide](../guides/estimation/observer-controller.md) -- composing observers with controllers
