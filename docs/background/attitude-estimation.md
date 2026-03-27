# Attitude Estimation

Attitude estimation determines the orientation of a rigid body in
three-dimensional space. Unlike Euclidean state estimation, orientations live
on the SO(3) manifold -- a curved space where standard vector addition does
not apply. Special formulations are needed to respect this geometry and avoid
singularities inherent in minimal rotation parameterisations
[1, Ch. 7, pp. 211--266].

Attitude estimation is fundamental to aerospace, robotics, and motion capture
systems, where IMU sensors (gyroscopes, accelerometers, magnetometers) provide
noisy measurements of angular velocity and reference field directions.

## Rotation Representations

Several representations exist for 3D rotations, each with trade-offs
[1, Sec. 7.1, pp. 213--222]:

- **Rotation matrix** ($3 \times 3$, 9 parameters, 6 constraints): direct,
  singularity-free, but redundant. The constraint $R^\top R = I$,
  $\det(R) = 1$ must be maintained numerically.
- **Quaternion** ($4$ parameters, $1$ constraint): compact, singularity-free,
  efficient composition via the Hamilton product. Unit quaternions form a
  double cover of SO(3): $q$ and $-q$ represent the same rotation.
- **Axis-angle** ($3$ parameters): the rotation axis $\hat{e}$ scaled by the
  rotation angle $\theta$, giving $\phi = \theta \hat{e}$. Minimal but the
  axis is undefined at zero rotation.
- **Euler angles** ($3$ parameters): intuitive but suffer from gimbal lock
  when the pitch angle reaches $\pm 90^\circ$.

For estimation, quaternions are the preferred representation due to their
compactness, lack of singularities, and efficient algebra
[2, Sec. 3.2, pp. 73--80].

## Quaternion Conventions

A unit quaternion $q = [q_w, q_x, q_y, q_z]^\top$ with
$\lVert q \rVert = 1$ represents a rotation. The Hamilton convention defines
the product such that quaternion multiplication corresponds to rotation
composition in the same order [1, Sec. 7.1.3, pp. 218--220].

The time evolution under body-frame angular velocity
$\omega = [\omega_x, \omega_y, \omega_z]^\top$ is
[1, Sec. 7.2, pp. 223--226]:

$$
\dot{q} = \frac{1}{2} q \otimes \begin{bmatrix} 0 \\ \omega \end{bmatrix}
$$

In discrete time with sampling period $T_s$, the first-order integration is:

$$
q_{k+1} = q_k \otimes \exp\!\left(\frac{T_s}{2} \omega_k\right)
$$

where $\otimes$ denotes the Hamilton quaternion product.

## Tangent Space and Lie Algebra

At any rotation on SO(3), the tangent space is a 3-dimensional vector space
isomorphic to $\mathbb{R}^3$ via the skew-symmetric matrix representation.
Small rotations are represented as 3-vectors in this tangent space using the
exponential and logarithmic maps [1, Sec. 7.3, pp. 227--234].

### Exponential Map (Rodrigues' Formula)

The exponential map takes a tangent vector $\phi \in \mathbb{R}^3$ to a
rotation [1, eq. (7.18), p. 229]:

$$
\exp([\phi]_\times) = I + \frac{\sin \theta}{\theta} [\phi]_\times + \frac{1 - \cos \theta}{\theta^2} [\phi]_\times^2
$$

where $\theta = \lVert \phi \rVert$ and $[\phi]_\times$ is the $3 \times 3$
skew-symmetric matrix. For quaternions:

$$
\exp(\phi) = \begin{bmatrix} \cos(\theta/2) \\ \frac{\sin(\theta/2)}{\theta} \phi \end{bmatrix}
$$

### Logarithmic Map

The inverse recovers the tangent vector from a quaternion
[1, eq. (7.22), p. 231]:

$$
\log(q) = \frac{2 \arccos(q_w)}{\lVert q_{xyz} \rVert} \, q_{xyz}
$$

with a sinc-expansion near $\theta = 0$ for numerical stability.

## Multiplicative EKF (MEKF)

The MEKF is the standard approach for quaternion-based attitude estimation.
It maintains the attitude estimate as a unit quaternion but performs Kalman
filtering on a 3-vector error state $\delta\theta$ in the tangent space. This
avoids the over-parameterisation of using a 4-dimensional quaternion state in
the filter [3, pp. 311--317].

### Error-State Definition

The attitude error is defined multiplicatively
[3, eq. (2), p. 312]:

$$
\delta q = \hat{q}^{-1} \otimes q_{\text{true}}
$$

For small errors, $\delta q \approx [1, \; \delta\theta/2]^\top$, giving the
3-vector error state $\delta\theta \in \mathbb{R}^3$.

### MEKF Prediction

Propagate the reference quaternion using the measured angular velocity
[1, Sec. 7.4, pp. 240--245]:

$$
\hat{q}_{k|k-1} = \hat{q}_{k-1|k-1} \otimes \exp\!\left(\frac{T_s}{2} \omega_k\right)
$$

Propagate the $3 \times 3$ error-state covariance:

$$
P_{k|k-1} = F_k \, P_{k-1|k-1} \, F_k^\top + G_k \, Q \, G_k^\top
$$

where $F_k \approx I - T_s [\omega_k]_\times$ is the discrete error-state
transition matrix and $G_k$ maps gyroscope noise to the error state.

### MEKF Update

Compute the Kalman gain for the 3-dimensional error state
[1, Sec. 7.4, pp. 245--248]:

$$
K_k = P_{k|k-1} \, H_k^\top (H_k \, P_{k|k-1} \, H_k^\top + R)^{-1}
$$

Apply the error-state correction:

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

The key insight is that the error state is unconstrained in $\mathbb{R}^3$,
so the standard Kalman filter applies without modification. The quaternion is
re-normalised after each correction to maintain the unit-norm constraint.

## Manifold UKF

The manifold UKF extends the unscented Kalman filter to SO(3) by performing
sigma-point operations in the tangent space [4, pp. 108--114]:

1. Generate sigma points in $\mathbb{R}^3$ from the error-state covariance
2. Map to the manifold: $q_i = \hat{q} \otimes \exp(\mathcal{X}_i)$
3. Propagate through nonlinear dynamics on the manifold
4. Compute statistics in the tangent space:
   $\delta_i = \log(\hat{q}^{-1} \otimes q_i')$
5. Recover mean and covariance from the tangent-space residuals

This preserves the geometry of SO(3) throughout the sigma-point transform
while avoiding the Jacobian computation required by the MEKF.

## Complementary Filter

For IMU-based attitude estimation, complementary filters offer a simpler and
computationally cheaper alternative to optimal filters. The Mahony filter
[5, pp. 1203--1218] fuses high-frequency gyroscope data (which drifts) with
low-frequency accelerometer/magnetometer data (which is noisy) using a
proportional-integral correction on SO(3):

$$
\omega_{\text{corr}} = \omega_{\text{gyro}} + K_p \, e_{\text{att}} + K_i \int e_{\text{att}} \, dt
$$

$$
\dot{q} = \frac{1}{2} q \otimes \begin{bmatrix} 0 \\ \omega_{\text{corr}} \end{bmatrix}
$$

The attitude error $e_{\text{att}}$ is computed from the cross product between
measured and predicted reference directions (gravity, magnetic north). The
gains $K_p$ and $K_i$ set the crossover frequency between gyroscope and
reference measurements.

Complementary filters are less optimal than the MEKF but require no covariance
propagation, no matrix inversions, and have deterministic constant-time
execution -- making them suitable for resource-constrained embedded systems.

## References

[1] T. D. Barfoot, "State Estimation for Robotics," 2nd ed., Cambridge
University Press, 2025.

[2] K. M. Lynch and F. C. Park, "Modern Robotics: Mechanics, Planning, and
Control," Cambridge University Press, 2017.

[3] F. L. Markley, "Attitude Error Representations for Kalman Filtering,"
Journal of Guidance, Control, and Dynamics, vol. 26, no. 2, pp. 311--317,
2003.

[4] S. Hauberg, F. Lauze, and K. S. Pedersen, "Unscented Kalman Filtering on
(Sub)Riemannian Manifolds," Journal of Mathematical Imaging and Vision, vol.
46, no. 1, pp. 103--120, 2013.

[5] R. Mahony, T. Hamel, and J.-M. Pflimlin, "Nonlinear Complementary
Filters on the Special Orthogonal Group," IEEE Transactions on Automatic
Control, vol. 53, no. 5, pp. 1203--1218, 2008.
