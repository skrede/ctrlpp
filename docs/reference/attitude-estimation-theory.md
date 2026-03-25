# Attitude Estimation Theory

Attitude estimation determines the orientation of a rigid body in
three-dimensional space. Unlike Euclidean state estimation, orientations live
on the SO(3) manifold -- a curved space where standard vector addition does
not apply. Special formulations are needed to handle the geometry correctly.

## Rotation Representations

Several representations exist for 3D rotations, each with trade-offs:

- **Rotation matrix** ($3 \times 3$, 9 parameters, 6 constraints): direct but
  redundant. Orthogonality must be maintained.
- **Quaternion** (4 parameters, 1 constraint): compact, singularity-free,
  efficient composition via Hamilton product. ctrlpp uses unit quaternions
  with Hamilton convention (w-first user-facing).
- **Axis-angle** (3 parameters): the rotation axis scaled by the rotation
  angle. Minimal, but has a singularity at zero rotation for the axis
  direction.
- **Euler angles** (3 parameters): intuitive but suffer from gimbal lock.
  Not used internally in ctrlpp.

## Quaternion Kinematics

A unit quaternion $q = [q_w, q_x, q_y, q_z]^\top$ with
$\lVert q \rVert = 1$ represents a rotation. The time evolution under angular
velocity $\omega = [\omega_x, \omega_y, \omega_z]^\top$ is:

$$
\dot{q} = \frac{1}{2} q \otimes \begin{bmatrix} 0 \\ \omega \end{bmatrix}
$$

In discrete time with sampling period $T_s$, the first-order integration is:

$$
q_{k+1} = q_k \otimes \exp\!\left(\frac{T_s}{2} \begin{bmatrix} 0 \\ \omega_k \end{bmatrix}\right)
$$

where $\otimes$ denotes the Hamilton quaternion product.

## The Tangent Space and Lie Algebra

At any point on SO(3), the tangent space is a 3-dimensional vector space
isomorphic to $\mathbb{R}^3$. Small rotations are represented as 3-vectors
in this tangent space using the exponential and logarithmic maps.

### Exponential Map (Rodrigues' Formula)

The exponential map takes a tangent vector $\phi \in \mathbb{R}^3$ to a
rotation matrix:

$$
\exp([\phi]_\times) = I + \frac{\sin \theta}{\theta} [\phi]_\times + \frac{1 - \cos \theta}{\theta^2} [\phi]_\times^2
$$

where $\theta = \lVert \phi \rVert$ and $[\phi]_\times$ is the
$3 \times 3$ skew-symmetric matrix of $\phi$:

$$
[\phi]_\times = \begin{bmatrix} 0 & -\phi_3 & \phi_2 \\ \phi_3 & 0 & -\phi_1 \\ -\phi_2 & \phi_1 & 0 \end{bmatrix}
$$

For quaternions, the exponential map is:

$$
\exp(\phi) = \begin{bmatrix} \cos(\theta/2) \\ \frac{\sin(\theta/2)}{\theta} \phi \end{bmatrix}
$$

### Logarithmic Map

The inverse operation recovers the tangent vector from a rotation:

$$
\log(q) = \frac{2 \arccos(q_w)}{\lVert q_{xyz} \rVert} \, q_{xyz}
$$

with a sinc-expansion near $\theta = 0$ for numerical stability.

## Multiplicative Error-State Formulation (MEKF)

The MEKF maintains the attitude estimate as a quaternion but performs Kalman
filtering on a 3-vector error state $\delta\theta$ in the tangent space.

### MEKF Prediction

Propagate the reference quaternion using angular velocity:

$$
\hat{q}_{k|k-1} = \hat{q}_{k-1|k-1} \otimes \exp\!\left(\frac{T_s}{2} \omega_k\right)
$$

Propagate the $3 \times 3$ error-state covariance:

$$
P_{k|k-1} = F_k \, P_{k-1|k-1} \, F_k^\top + G_k \, Q \, G_k^\top
$$

where $F_k \approx I - T_s [\omega_k]_\times$ is the discrete error-state
transition matrix.

### MEKF Update

Compute the Kalman gain for the 3D error state:

$$
K_k = P_{k|k-1} \, H_k^\top (H_k \, P_{k|k-1} \, H_k^\top + R)^{-1}
$$

Apply the error-state correction to the quaternion:

$$
\delta\theta_k = K_k \, (z_k - h(\hat{q}_{k|k-1}))
$$

$$
\hat{q}_{k|k} = \hat{q}_{k|k-1} \otimes \exp(\delta\theta_k / 2)
$$

Reset the error state to zero and update the covariance:

$$
P_{k|k} = (I - K_k H_k) \, P_{k|k-1}
$$

This avoids the quaternion unit-norm constraint inside the filter -- the
error state is unconstrained in $\mathbb{R}^3$. The quaternion is always
re-normalised after the correction.

## Manifold UKF

The manifold UKF extends the unscented Kalman filter to SO(3) by:

1. Generating sigma points in the tangent space $\mathbb{R}^3$
2. Mapping them to the manifold via the exponential map:
   $q_i = \hat{q} \otimes \exp(\mathcal{X}_i)$
3. Propagating through the dynamics on the manifold
4. Computing statistics back in the tangent space via the logarithmic map:
   $\delta_i = \log(\hat{q}^{-1} \otimes q_i')$

This preserves the geometry of SO(3) throughout the sigma-point transform.

## Complementary Filtering

For IMU-based attitude estimation, complementary filters fuse gyroscope
(high-frequency, drifts) with accelerometer/magnetometer (low-frequency,
noisy). The Mahony filter uses a proportional-integral correction on SO(3):

$$
\omega_{\text{corr}} = \omega_{\text{gyro}} + K_p \, e_{\text{att}} + K_i \int e_{\text{att}} \, dt
$$

$$
\dot{q} = \frac{1}{2} q \otimes \begin{bmatrix} 0 \\ \omega_{\text{corr}} \end{bmatrix}
$$

The attitude error $e_{\text{att}}$ is computed from the cross product between
measured and predicted gravity/magnetic field directions. Complementary
filters are simpler and cheaper than MEKF but less optimal.

## References

- Sola, J., Deray, J., and Atchuthan, D. (2018). "A Micro Lie Theory for
  State Estimation in Robotics." arXiv:1812.01537. [`sola2018`]
  Concise and practical introduction to Lie groups for robotics, covering
  SO(3), SE(3), and their use in estimation. Primary reference for ctrlpp's
  SO(3) implementation.

- Markley, F. L. (2003). "Attitude Error Representations for Kalman
  Filtering." *J. Guidance, Control, and Dynamics*, 26(2):311--317.
  [`markley2003`]
  Establishes the multiplicative error-state formulation for quaternion-based
  attitude estimation.

- Mahony, R., Hamel, T., and Pflimlin, J.-M. (2008). "Nonlinear
  Complementary Filters on the Special Orthogonal Group." *IEEE Trans.
  Automatic Control*, 53(5):1203--1218. [`mahony2008`]
  The foundational paper for complementary attitude filtering on SO(3).

- Barfoot, T. D. (2017). *State Estimation for Robotics.* Cambridge
  University Press. [`barfoot2017`]
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
