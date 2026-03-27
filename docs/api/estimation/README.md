# Estimation

State estimation and observer types for inferring hidden system states from noisy
measurements. The family spans linear observers (Kalman, Luenberger), nonlinear
filters (EKF, UKF, particle filter), and attitude-aware estimators that operate
directly on SO(3) (MEKF, manifold UKF, complementary filter).

## Types

### Linear

- [kalman](kalman.md) -- Linear Kalman filter (predict-update cycle)
- [luenberger](luenberger.md) -- Luenberger state observer (fixed gain)

### Nonlinear

- [ekf](ekf.md) -- Extended Kalman filter (analytical or numerical Jacobians)
- [ukf](ukf.md) -- Unscented Kalman filter (sigma point strategies)
- [particle_filter](particle-filter.md) -- Bootstrap SIR particle filter (ESS-adaptive resampling)

### Attitude-Aware

- [mekf](mekf.md) -- Multiplicative extended Kalman filter for SO(3) attitude
- [manifold_ukf](manifold-ukf.md) -- Manifold unscented Kalman filter for SO(3) attitude
- [complementary_filter](complementary-filter.md) -- Mahony complementary filter (IMU/MARG)

### Interface

- [observer_policy](observer-policy.md) -- Observer concept for controller composition

## When to use

Pick **kalman** or **luenberger** for linear systems with known dynamics.

Pick **EKF** when you have a nonlinear model and can provide (or auto-generate)
Jacobians. Pick **UKF** when Jacobians are unavailable or the nonlinearity is
severe -- sigma points handle it without linearisation.

Pick **particle filter** for highly nonlinear or multimodal distributions where
Gaussian assumptions break down.

Pick **MEKF** or **manifold UKF** for attitude estimation where states live on
SO(3). Pick **complementary filter** for lightweight IMU/MARG sensor fusion.

## Theory

For the mathematical background, see
[Kalman Theory](../../background/kalman.md),
[EKF/UKF Theory](../../background/ekf-ukf.md),
[Particle Filter Theory](../../background/particle-filter.md), and
[Attitude Estimation Theory](../../background/attitude-estimation.md).
