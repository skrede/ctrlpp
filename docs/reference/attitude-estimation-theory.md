# Attitude Estimation Theory

Attitude estimation determines the orientation of a rigid body in
three-dimensional space. Unlike Euclidean state estimation, orientations live
on the SO(3) manifold -- a curved space where standard vector addition does
not apply. Special formulations are needed to handle the geometry correctly.

## Key Concepts

### Rotation Representations

Several representations exist for 3D rotations, each with trade-offs:

- **Rotation matrix** (3x3, 9 parameters, 6 constraints): direct but
  redundant. Orthogonality must be maintained.
- **Quaternion** (4 parameters, 1 constraint): compact, singularity-free,
  efficient composition via Hamilton product. ctrlpp uses unit quaternions
  with Hamilton convention (w-first user-facing).
- **Axis-angle** (3 parameters): the rotation axis scaled by the rotation
  angle. Minimal, but has a singularity at zero rotation for the axis
  direction.
- **Euler angles** (3 parameters): intuitive but suffer from gimbal lock.
  Not used internally in ctrlpp.

### The Tangent Space

At any point on SO(3), the tangent space is a 3-dimensional vector space
isomorphic to R^3. Small rotations are represented as 3-vectors in this
tangent space using the exponential map:

```
R = exp([phi]_x)     (tangent vector phi -> rotation R)
phi = log(R)          (rotation R -> tangent vector phi)
```

where [phi]_x is the skew-symmetric matrix of phi. This is the Lie algebra
so(3).

### Multiplicative Error-State Formulation (MEKF)

The MEKF maintains the attitude estimate as a quaternion but performs Kalman
filtering on a 3-vector error state in the tangent space:

1. **Predict**: propagate the quaternion using angular velocity, propagate
   the 3x3 error-state covariance
2. **Update**: compute the Kalman gain for the 3D error state, compute the
   error-state correction, apply it to the quaternion via the exponential
   map, reset the error state to zero

This avoids the quaternion unit-norm constraint inside the filter -- the
error state is unconstrained in R^3. The quaternion is always re-normalised
after the correction.

### Manifold UKF

The manifold UKF extends the unscented Kalman filter to manifolds by:

1. Generating sigma points in the tangent space
2. Mapping them to the manifold via the exponential map
3. Propagating through the dynamics on the manifold
4. Computing statistics back in the tangent space via the logarithmic map

This preserves the geometry of SO(3) throughout the sigma-point transform.

### Complementary Filtering

For IMU-based attitude estimation, complementary filters fuse gyroscope
(high-frequency, drifts) with accelerometer/magnetometer (low-frequency,
noisy). The Mahony filter uses a proportional-integral correction on SO(3):

```
omega_corrected = omega_gyro + Kp * error + Ki * integral(error)
q_dot = 0.5 * q * [0, omega_corrected]
```

The error is computed from the cross product between measured and predicted
gravity/magnetic field directions. Complementary filters are simpler and
cheaper than MEKF but less optimal.

## References

- **Sola, J., Deray, J., and Atchuthan, D.** "A Micro Lie Theory for State
  Estimation in Robotics." arXiv:1812.01537, 2018. Published in *Int. J.
  Robotics Research*, 2021.
  Concise and practical introduction to Lie groups for robotics, covering
  SO(3), SE(3), and their use in estimation. Primary reference for ctrlpp's
  SO(3) implementation.

- **Markley, F. L.** "Attitude Error Representations for Kalman Filtering."
  *J. Guidance, Control, and Dynamics*, 26(2):311--317, 2003.
  DOI: 10.2514/2.5048.
  Establishes the multiplicative error-state formulation for quaternion-based
  attitude estimation.

- **Mahony, R., Hamel, T., and Pflimlin, J.-M.** "Nonlinear Complementary
  Filters on the Special Orthogonal Group." *IEEE Trans. Automatic Control*,
  53(5):1203--1218, 2008. DOI: 10.1109/TAC.2008.923738.
  The foundational paper for complementary attitude filtering on SO(3).

- **Barfoot, T. D.** *State Estimation for Robotics.* Cambridge University
  Press, 2017. ISBN 978-1-107-15939-2.
  Comprehensive textbook covering Lie group estimation, including MEKF and
  manifold sigma-point methods.

## Related API Pages

- [so3](../lie/so3.md) -- quaternion SO(3) utilities (exp, log, compose,
  skew)
- [mekf](../estimation/mekf.md) -- multiplicative extended Kalman filter
- [manifold_ukf](../estimation/manifold-ukf.md) -- manifold unscented
  Kalman filter
- [complementary_filter](../estimation/complementary-filter.md) -- Mahony
  complementary filter
