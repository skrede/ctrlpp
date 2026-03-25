# Particle Filter Theory

Particle filters (sequential Monte Carlo methods) represent the posterior
distribution of the state using a set of weighted samples (particles). Unlike
the Kalman family, particle filters make no assumptions about linearity or
Gaussianity, making them applicable to arbitrary distributions and strongly
nonlinear systems.

## Key Concepts

### Importance Sampling

The posterior distribution p(x[k] | z[1:k]) is approximated by N weighted
particles:

```
p(x[k] | z[1:k]) ~ sum_i  w_i * delta(x - x_i)
```

Each particle represents a hypothesis about the true state. Particles with
higher weights are more consistent with the observations.

### Bootstrap Filter (SIR)

The Sequential Importance Resampling (SIR) or bootstrap filter is the
simplest particle filter variant:

1. **Propagate**: draw each particle through the dynamics model with process
   noise: x_i[k] ~ p(x[k] | x_i[k-1], u[k-1])
2. **Weight**: compute the likelihood of the measurement given each particle:
   w_i = p(z[k] | x_i[k])
3. **Normalise**: scale weights so they sum to one
4. **Resample**: draw N particles with replacement, proportional to weights

This is the variant implemented in ctrlpp.

### Resampling

After several steps, most particles accumulate negligible weight while a few
particles carry all the mass. This *weight degeneracy* reduces the effective
sample size. Resampling duplicates high-weight particles and discards
low-weight ones, restoring diversity.

The effective sample size (ESS) measures degeneracy:

```
ESS = 1 / sum(w_i^2)
```

When ESS drops below a threshold (typically N/2), resampling is triggered.
ctrlpp uses ESS-adaptive resampling: resample only when needed.

### Particle Diversity and Roughening

After resampling, many particles may be identical copies. *Roughening* adds
small random perturbations to the resampled particles, preventing sample
impoverishment and maintaining diversity.

### Convergence

As the number of particles N increases, the particle approximation converges
to the true posterior. In practice, the required N depends on the state
dimension -- particle filters become expensive in high dimensions (the "curse
of dimensionality"). They are most effective for low-dimensional state spaces
(typically fewer than 6--10 dimensions).

## References

- **Gordon, N. J., Salmond, D. J., and Smith, A. F. M.** "Novel Approach to
  Nonlinear/Non-Gaussian Bayesian State Estimation." *IEE Proceedings F*,
  140(2):107--113, 1993. DOI: 10.1049/ip-f-2.1993.0015.
  The foundational paper introducing the bootstrap particle filter.

- **Arulampalam, M. S., Maskell, S., Gordon, N., and Clapp, T.** "A Tutorial
  on Particle Filters for Online Nonlinear/Non-Gaussian Bayesian Tracking."
  *IEEE Transactions on Signal Processing*, 50(2):174--188, 2002.
  DOI: 10.1109/78.978374.
  Comprehensive tutorial covering SIR, auxiliary particle filters, and
  practical implementation considerations.

## Related API Pages

- [particle_filter](../estimation/particle-filter.md) -- bootstrap SIR
  particle filter with compile-time particle count and ESS-adaptive
  resampling
- [EKF/UKF Theory](ekf-ukf-theory.md) -- alternative nonlinear estimators
- [Observer-Controller Guide](../guides/estimation/observer-controller.md) -- composing observers with controllers
