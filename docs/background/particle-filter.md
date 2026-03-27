# Particle Filter

Particle filters, also known as sequential Monte Carlo (SMC) methods,
represent the posterior distribution of the state using a set of weighted
random samples called particles. Unlike the Kalman family, particle filters
make no assumptions about linearity or Gaussianity, making them applicable to
arbitrary distributions and strongly nonlinear systems
[1, Ch. 15, pp. 461--488].

The bootstrap particle filter was introduced by Gordon, Salmond, and Smith
in 1993 [2] and has since become the standard approach for low-dimensional
nonlinear/non-Gaussian estimation problems.

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
[1, Sec. 15.2, pp. 465--468]:

$$
p(x_k \mid z_{1:k}) \approx \sum_{i=1}^{N_p} w_k^{(i)} \, \delta(x_k - x_k^{(i)})
$$

where $x_k^{(i)}$ is the $i$-th particle, $w_k^{(i)}$ is its normalised
weight, and $\delta$ is the Dirac delta. Each particle represents a hypothesis
about the true state; particles with higher weights are more consistent with
the observations.

The concept derives from importance sampling: samples are drawn from a
proposal distribution $q(x_k \mid x_{k-1}, z_k)$ and their weights are
adjusted to account for the discrepancy between the proposal and the true
posterior [3, Sec. 2, pp. 176--178].

## Bootstrap Filter (SIR)

The Sequential Importance Resampling (SIR) or bootstrap filter is the simplest
and most widely used variant. It uses the prior
$p(x_k \mid x_{k-1}^{(i)}, u_{k-1})$ as the proposal distribution
[2, Sec. 3, pp. 109--111].

### Algorithm

**1. Propagate** -- draw each particle through the dynamics model with process
noise [1, Sec. 15.3, pp. 469--472]:

$$
x_k^{(i)} \sim p(x_k \mid x_{k-1}^{(i)}, u_{k-1})
$$

In practice: $x_k^{(i)} = f(x_{k-1}^{(i)}, u_{k-1}) + w_k^{(i)}$ where
$w_k^{(i)}$ is drawn from the process noise distribution.

**2. Weight** -- compute the likelihood of the measurement given each particle
[1, Sec. 15.3, pp. 472--474]:

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

**4. Resample** -- when triggered by the ESS criterion, draw $N_p$ particles
with replacement, proportional to weights.

### State Estimate

The weighted mean provides the minimum mean square error (MMSE) estimate
[1, Sec. 15.3, p. 474]:

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
the effective sample size [3, Sec. 3.1, pp. 179--180]:

$$
N_{\text{eff}} = \frac{1}{\sum_{i=1}^{N_p} \left(w_k^{(i)}\right)^2}
$$

When $N_{\text{eff}}$ drops below a threshold (typically $N_p / 2$),
resampling is triggered.

### Systematic Resampling

Systematic resampling generates a single random offset
$u_0 \sim \mathcal{U}(0, 1/N_p)$ and takes samples at
$u_0 + j/N_p$ for $j = 0, \ldots, N_p - 1$ from the cumulative weight
distribution [3, Sec. 3.2, pp. 181--182]. This has lower variance than
multinomial resampling and is $O(N_p)$.

## Roughening

After resampling, many particles may be identical copies, leading to sample
impoverishment. Roughening adds small random perturbations to maintain
diversity [1, Sec. 15.4, pp. 478--480]:

$$
x_k^{(i)} \leftarrow x_k^{(i)} + \epsilon^{(i)}, \qquad \epsilon^{(i)} \sim \mathcal{N}(0, \sigma_r^2 I)
$$

where $\sigma_r$ is proportional to the spread of the particle cloud. This
prevents the particle set from collapsing to a single point.

## Convergence and Curse of Dimensionality

As $N_p \to \infty$, the particle approximation converges to the true
posterior [1, Sec. 15.5, pp. 481--484]. The convergence rate is
$O(1/\sqrt{N_p})$ independent of dimension. However, the number of particles
required for a given accuracy grows exponentially with state dimension -- the
"curse of dimensionality" [3, Sec. 4, pp. 183--185].

In practice, particle filters are most effective for low-dimensional state
spaces (typically fewer than 6--10 dimensions). For higher dimensions,
Rao-Blackwellisation can be used to marginalise out linear sub-states
analytically [1, Sec. 15.6, pp. 485--488].

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
