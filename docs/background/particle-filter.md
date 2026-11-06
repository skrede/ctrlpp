# Particle Filter

Particle filters, also known as sequential Monte Carlo (SMC) methods,
represent the posterior distribution of the state using a set of weighted
random samples called particles. Unlike the Kalman family, particle filters
make no assumptions about linearity or Gaussianity, making them applicable to
arbitrary distributions and strongly nonlinear systems
[1, Ch. 15, pp. 461--483].

The bootstrap particle filter was introduced by Gordon, Salmond, and Smith
in 1993 [2]. ctrlpp's
[particle_filter](../api/estimation/particle-filter.md) implements it with
narrower noise models than the general setting below; the section
[What ctrlpp implements](#what-ctrlpp-implements) states them.

## Problem Setting

Consider a discrete nonlinear system with possibly non-Gaussian noise:

$$
x_{k+1} = f(x_k, u_k, w_k)
$$

$$
z_k = h(x_k, v_k)
$$

where the noise $w_k$ and $v_k$ may have arbitrary distributions. The
objective is to approximate the posterior distribution
$p(x_k \mid z_{1:k})$ at each time step [1, Sec. 15.1, pp. 462--464].

## Importance Sampling

The posterior distribution is approximated by $N_p$ weighted particles
[3, Sec. V-A, Eq. (49), p. 178]:

$$
p(x_k \mid z_{1:k}) \approx \sum_{i=1}^{N_p} w_k^{(i)} \, \delta(x_k - x_k^{(i)})
$$

where $x_k^{(i)}$ is the $i$-th particle, $w_k^{(i)}$ is its normalized
weight, and $\delta$ is the Dirac delta. Each particle represents a hypothesis
about the true state; particles with higher weights are more consistent with
the observations.

The concept derives from importance sampling: samples are drawn from a
proposal distribution $q(x_k \mid x_{k-1}, z_k)$ and their weights are
adjusted to account for the discrepancy between the proposal and the true
posterior [3, Sec. V-A, Eqs. (41)--(48), pp. 177--178].

## Bootstrap Filter (SIR)

The Sequential Importance Resampling (SIR) or bootstrap filter is the simplest
and most widely used variant. It uses the prior
$p(x_k \mid x_{k-1}^{(i)}, u_{k-1})$ as the proposal distribution
[2, Sec. 3.1, p. 108].

### Algorithm

**1. Propagate**<br/>draw each particle through the dynamics model with process
noise [1, Sec. 15.2, Eq. (15.17), p. 466]:

$$
x_k^{(i)} \sim p(x_k \mid x_{k-1}^{(i)}, u_{k-1})
$$

In practice: $x_k^{(i)} = f(x_{k-1}^{(i)}, u_{k-1}) + w_k^{(i)}$ where
$w_k^{(i)}$ is drawn from the process noise distribution.

**2. Weight**<br/>compute the likelihood of the measurement given each particle
[1, Sec. 15.2, Eq. (15.18), p. 466]:

$$
\tilde{w}_k^{(i)} = p(z_k \mid x_k^{(i)})
$$

When resampling is not done every step, as with the ESS criterion below, the
previous weight is carried:
$\tilde{w}_k^{(i)} = w_{k-1}^{(i)}\, p(z_k \mid x_k^{(i)})$ [3, Eq. (63), p. 179]. ctrlpp
does this.

For Gaussian measurement noise $v_k \sim \mathcal{N}(0, R)$:

$$
\tilde{w}_k^{(i)} \propto \exp\!\left(-\frac{1}{2} (z_k - h(x_k^{(i)}))^\top R^{-1} (z_k - h(x_k^{(i)}))\right)
$$

**3. Normalize**<br/>scale weights to sum to one:

$$
w_k^{(i)} = \frac{\tilde{w}_k^{(i)}}{\sum_{j=1}^{N_p} \tilde{w}_k^{(j)}}
$$

**4. Resample**<br/>when triggered by the ESS criterion, draw $N_p$ particles
with replacement, proportional to weights, and reset every weight to $1/N_p$.

### State Estimate

The weighted mean provides the minimum mean square error (MMSE) estimate
[1, Sec. 15.2, Eq. (15.21), p. 468]:

$$
\hat{x}_k = \sum_{i=1}^{N_p} w_k^{(i)} \, x_k^{(i)}
$$

The weighted covariance quantifies estimation uncertainty:

$$
P_k \approx \sum_{i=1}^{N_p} w_k^{(i)} (x_k^{(i)} - \hat{x}_k)(x_k^{(i)} - \hat{x}_k)^\top
$$

## Effective Sample Size and Resampling

After several steps without resampling, most particles accumulate negligible
weight while a few carry all the mass. This weight degeneracy is measured by
the effective sample size [3, Sec. V-A, Eq. (51), p. 179]:

$$
N_{\text{eff}} = \frac{1}{\sum_{i=1}^{N_p} \left(w_k^{(i)}\right)^2}
$$

When $N_{\text{eff}}$ drops below a threshold $N_T$, resampling is triggered
[3, Sec. V-A, Alg. 3, p. 180]. ctrlpp's default is $N_T = N_p / 2$, set through
`ess_threshold`.

### Systematic Resampling

Systematic resampling generates a single random offset
$u_0 \sim \mathcal{U}(0, 1/N_p)$ and takes samples at
$u_0 + j/N_p$ for $j = 0, \ldots, N_p - 1$ from the cumulative weight
distribution [3, Sec. V-A, Alg. 2, p. 180]. This has lower variance than
multinomial resampling and is $O(N_p)$. It is ctrlpp's default; multinomial
resampling, $N_p$ independent uniform draws, is also provided.

## Roughening

After resampling, many particles may be identical copies, leading to sample
impoverishment. Roughening adds small random perturbations to maintain
diversity [1, Sec. 15.3.1.1, Eq. (15.26), pp. 470--471]:

$$
x_k^{(i)} \leftarrow x_k^{(i)} + \epsilon^{(i)}, \qquad \epsilon^{(i)} \sim \mathcal{N}\!\left(0, \operatorname{diag}(\sigma_1^2, \ldots, \sigma_{n_x}^2)\right)
$$

where each $\sigma_d$ is proportional to the spread of the particle cloud in
coordinate $d$. This prevents the particle set from collapsing to a single
point. ctrlpp uses $\sigma_d = K E_d N_p^{-1/n_x}$ [2, Sec. 4.2, p. 112], with $E_d$ the range
(maximum minus minimum) of coordinate $d$ over the resampled set and $K$ the
configured `roughening_scale`; a coordinate with zero range falls back to
$\sigma_d = K \sqrt{Q_{dd}}$. Roughening is applied only after a resampling.

## Convergence and Curse of Dimensionality

As $N_p \to \infty$, the particle approximation converges to the true
posterior [1, Sec. 15.2, p. 467]. The Monte Carlo error is
$O(1/\sqrt{N_p})$ independent of dimension, but only in theory and not
necessarily in practice: the number of particles required for a given
accuracy still grows with the state dimension &mdash; the "curse of
dimensionality" [4, Sec. 11.1].

In practice, particle filters are most effective for low-dimensional state
spaces (typically fewer than 6--10 dimensions). For higher dimensions,
Rao-Blackwellization can be used to marginalize out linear sub-states
analytically [4, Sec. 11.7]; ctrlpp does not implement it.

## What ctrlpp implements

- The process model is $x_k = f(x_{k-1}, u_{k-1}) + w_k$ with
  $w_k \sim \mathcal{N}(0, Q)$: the dynamics callable takes no noise input,
  and the noise is additive Gaussian with a discrete covariance $Q$ drawn
  once per particle per step. The general $f(x_k, u_k, w_k)$ above is not
  supported.
- The default likelihood is the Gaussian above with covariance $R$, including
  its normalizer, on the plain difference $z_k - h(x_k^{(i)})$. Any other
  measurement density, including one on a wrapped domain such as a bearing, is
  supplied as a likelihood policy returning
  $\log p(z_k \mid x_k^{(i)})$; see the
  [likelihood policy](../api/estimation/particle-filter.md).
- The prior is $\mathcal{N}(x_0, P_0)$.
- Weights are kept as logarithms (normalized with log-sum-exp) or as linear
  values. When no particle explains a measurement, the weights are reset to
  uniform and the filter reports it.
- The estimate is the weighted mean above or the highest-weight particle; no
  mode estimate is formed.

## References

[1] D. Simon, "Optimal State Estimation: Kalman, H-Infinity, and Nonlinear
Approaches," Wiley, 2006.

[2] N. J. Gordon, D. J. Salmond, and A. F. M. Smith, "Novel Approach to
Nonlinear/Non-Gaussian Bayesian State Estimation," IEE Proceedings F, vol.
140, no. 2, pp. 107--113, 1993.

[3] M. S. Arulampalam, S. Maskell, N. Gordon, and T. Clapp, "A Tutorial on
Particle Filters for Online Nonlinear/Non-Gaussian Bayesian Tracking," IEEE
Transactions on Signal Processing, vol. 50, no. 2, pp. 174--188, 2002.

[4] S. Sarkka and L. Svensson, "Bayesian Filtering and Smoothing," 2nd ed.,
Cambridge University Press, 2023.
