# Attitude Estimation

Attitude estimation determines the orientation of a rigid body in
three-dimensional space. Unlike Euclidean state estimation, orientations live
on the SO(3) manifold &mdash; a curved space where standard vector addition does
not apply. Special formulations are needed to respect this geometry and avoid
singularities inherent in minimal rotation parametrizations
[1, Secs. 7.2.3 and 8.1.1].

Attitude estimation is fundamental to aerospace, robotics, and motion capture
systems, where IMU sensors (gyroscopes, accelerometers, magnetometers) provide
noisy measurements of angular velocity and reference field directions.

## Rotation Representations

Several representations exist for 3D rotations, each with trade-offs
[1, Secs. 7.2.1 and 7.2.3]:

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

For estimation, unit quaternions avoid the singularities of every
three-parameter representation and form the rotation matrix without
trigonometric functions, at the cost of a fourth parameter and a unit-norm
constraint the estimator must enforce [1, Sec. 7.2.3], [2, App. B.3].

## Quaternion Conventions

A unit quaternion $q = [q_w, q_x, q_y, q_z]^\top$ with
$\lVert q \rVert = 1$ represents a rotation. The Hamilton convention defines
the product such that quaternion multiplication corresponds to rotation-matrix
multiplication in the same order, $R(p \otimes q) = R(p)\,R(q)$
[2, App. B.3]. ctrlpp uses this convention with the w-first
order shown.

The time evolution under body-frame angular velocity
$\omega = [\omega_x, \omega_y, \omega_z]^\top$ is
[5, App. B, Eq. (45)]:

$$
\dot{q} = \frac{1}{2} q \otimes \begin{bmatrix} 0 \\ \omega \end{bmatrix}
$$

Holding $\omega_k$ constant over a sampling period $T_s$, the exact solution
of this equation over the step is:

$$
q_{k+1} = q_k \otimes \exp\!\left(T_s\, \omega_k\right)
$$

where $\otimes$ denotes the Hamilton quaternion product and $\exp$ is the
quaternion exponential below, which takes the full rotation vector
$T_s \omega_k$ and halves the angle itself. This is the integration every
ctrlpp attitude filter performs, followed by a renormalization.

## Tangent Space and Lie Algebra

At any rotation on SO(3), the tangent space is a 3-dimensional vector space
isomorphic to $\mathbb{R}^3$ via the skew-symmetric matrix representation.
Small rotations are represented as 3-vectors in this tangent space using the
exponential and logarithmic maps [1, Secs. 8.1.2 and 8.1.3].

### Exponential Map (Rodrigues' Formula)

The exponential map takes a tangent vector $\phi \in \mathbb{R}^3$ to a
rotation [1, Sec. 8.1.3, Eq. (8.23)]:

$$
\exp([\phi]_\times) = I + \frac{\sin \theta}{\theta} [\phi]_\times + \frac{1 - \cos \theta}{\theta^2} [\phi]_\times^2
$$

where $\theta = \lVert \phi \rVert$ and $[\phi]_\times$ is the $3 \times 3$
skew-symmetric matrix. For quaternions [6, App. B, Eq. (132)]:

$$
\exp(\phi) = \begin{bmatrix} \cos(\theta/2) \\ \frac{\sin(\theta/2)}{\theta} \phi \end{bmatrix}
$$

### Logarithmic Map

The inverse recovers the tangent vector from a quaternion
[6, App. B, Eq. (133)]. Because $q$ and $-q$ are the same rotation, ctrlpp
first flips $q$ to the hemisphere $q_w \ge 0$, then computes

$$
\log(q) = \frac{2\, \operatorname{atan2}(\lVert q_{xyz} \rVert, q_w)}{\lVert q_{xyz} \rVert} \, q_{xyz}
$$

so the result is unique and its norm, the rotation angle, lies in
$[0, \pi]$. The `atan2` form equals $2 \arccos(q_w)$ for a unit quaternion
but keeps its accuracy near $\theta = 0$, where $\arccos$ loses it. Below
$\lVert q_{xyz} \rVert = 10^{-7}$ the limit $\log(q) = 2\, q_{xyz}$ is used.
The filters align each quaternion with the reference before taking a
difference, $\log(\hat{q}^{-1} \otimes q_i)$ with $q_i$ negated when
$\hat{q} \cdot q_i < 0$, so a sigma point and its antipode give the same
residual.

## Multiplicative EKF (MEKF)

The MEKF is the standard approach for quaternion-based attitude estimation.
It maintains the attitude estimate as a unit quaternion but performs Kalman
filtering on a 3-vector error state $\delta\theta$ in the tangent space. This
avoids the over-parametrization of using a 4-dimensional quaternion state in
the filter [3, pp. 311--317]. ctrlpp's MEKF also estimates a bias
$\hat{b} \in \mathbb{R}^{N_B}$, $N_B \ge 3$, whose first three components
are the gyroscope bias, so its error state is
$[\delta\theta; \delta b] \in \mathbb{R}^{3 + N_B}$.

### Error-State Definition

The attitude error is defined multiplicatively
[3, eq. (2), p. 312]:

$$
\delta q = \hat{q}^{-1} \otimes q_{\text{true}}
$$

that is, $q_{\text{true}} = \hat{q} \otimes \exp(\delta\theta)$, a
body-frame (right) error. For small errors,
$\delta q \approx [1, \; \delta\theta/2]^\top$, giving the 3-vector error
state $\delta\theta \in \mathbb{R}^3$.

### MEKF Prediction

Propagate the reference quaternion using the bias-corrected angular velocity
$\bar{\omega}_k = \omega_k - \hat{b}_{1:3}$
[3]:

$$
\hat{q}_{k|k-1} = \hat{q}_{k-1|k-1} \otimes \exp\!\left(T_s\, \bar{\omega}_k\right)
$$

The bias is held. Propagate the $(3 + N_B) \times (3 + N_B)$ error-state
covariance:

$$
P_{k|k-1} = F_k \, P_{k-1|k-1} \, F_k^\top + Q,
\qquad
F_k = \begin{bmatrix} C_k^\top & -T_s \left[ I_3 \;\; 0 \right] \\ 0 & I_{N_B} \end{bmatrix}
$$

where $C_k$ is the rotation matrix of $\exp(T_s \bar{\omega}_k)$, so
$C_k^\top \approx I - T_s [\bar{\omega}_k]_\times$ to first order. The
attitude-bias block is the first-order term $-T_s I$; the exact discrete
coefficient carries the right Jacobian of $T_s \bar{\omega}_k$, which ctrlpp
does not form. $Q$ is a discrete covariance added once per step as given:
ctrlpp has no noise-input matrix and does not scale $Q$ by $T_s$.

### MEKF Update

Compute the Kalman gain for the error state
[3], with $H_k$ the Jacobian of $h$ with respect to
$[\delta\theta; \delta b]$ at the nominal state:

$$
K_k = P_{k|k-1} \, H_k^\top (H_k \, P_{k|k-1} \, H_k^\top + R)^{-1}
$$

Apply the error-state correction:

$$
\begin{bmatrix} \delta\theta_k \\ \delta b_k \end{bmatrix} = K_k \, (z_k - h(\hat{q}_{k|k-1}, \hat{b}_{k|k-1}))
$$

$$
\hat{q}_{k|k} = \hat{q}_{k|k-1} \otimes \exp(\delta\theta_k),
\qquad
\hat{b}_{k|k} = \hat{b}_{k|k-1} + \delta b_k
$$

The quaternion is renormalized after the correction. Update the covariance
in Joseph form:

$$
P^{+} = (I - K_k H_k) \, P_{k|k-1} \, (I - K_k H_k)^\top + K_k R K_k^\top
$$

Folding $\delta\theta_k$ into the reference moves the frame the error state
is expressed in, so after the error state is reset to zero the covariance is
carried into the new frame:

$$
P_{k|k} = G_k \, P^{+} \, G_k^\top,
\qquad
G_k = \begin{bmatrix} I - \tfrac{1}{2} [\delta\theta_k]_\times & 0 \\ 0 & I \end{bmatrix}
$$

The error state is unconstrained in $\mathbb{R}^{3 + N_B}$, so the Kalman
update itself is the standard one; the multiplicative correction and this
reset are what the manifold adds.

## Manifold UKF

The manifold UKF extends the unscented Kalman filter to SO(3) by performing
sigma-point operations in the tangent space [4, Sec. 3]:

1. Generate sigma points $\mathcal{X}_i$ around zero in $\mathbb{R}^3$ from
   the tangent-space covariance
2. Map to the manifold: $q_i = \hat{q} \otimes \exp(\mathcal{X}_i)$
3. Propagate through nonlinear dynamics on the manifold, $q_i' = f(q_i, \omega)$
4. Find the weighted geodesic (Karcher) mean $\bar{q}$ by iterating
   $\bar{q} \leftarrow \bar{q} \otimes \exp\!\left(\sum_i W_i^m \log(\bar{q}^{-1} \otimes q_i')\right)$
   from $q_0'$ until the step's norm falls below a tolerance or an iteration
   budget runs out
5. Recover the covariance from the tangent-space residuals
   $\delta_i = \log(\bar{q}^{-1} \otimes q_i')$ about that mean,
   $P = \sum_i W_i^c \delta_i \delta_i^\top + Q$

The update draws fresh sigma points around the predicted attitude, forms
$\hat{z} = \sum_i W_i^m h(q_i)$, the innovation covariance $S$ (plus $R$)
and the cross-covariance from the residuals $\log(\hat{q}^{-1} \otimes q_i)$,
corrects $\hat{q} \leftarrow \hat{q} \otimes \exp(K (z - \hat{z}))$ with
$K = P_{xz} S^{-1}$, and carries $P - K S K^\top$ into the corrected frame with
the same reset Jacobian as the MEKF. ctrlpp's manifold UKF estimates the
attitude alone, with no bias state.

This preserves the geometry of SO(3) throughout the sigma-point transform
while needing no measurement Jacobian.

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

ctrlpp carries the integral as a bias estimate $\hat{b} = -K_i \int e\,dt$
and takes each step as:

$$
\hat{b}_k = \hat{b}_{k-1} - K_i \, e_k \, T_s,
\qquad
q_k = q_{k-1} \otimes \exp\!\left(T_s \left(\omega_{\text{gyro}} - \hat{b}_k + K_p \, e_k\right)\right)
$$

with $e_k = \hat{a} \times \hat{R}^\top e_z$ for the accelerometer direction
$\hat{a}$, plus, with a magnetometer direction $\hat{m}$,
$\hat{m} \times \hat{R}^\top b$ where $h = \hat{R}\hat{m}$ and
$b = [\sqrt{h_x^2 + h_y^2},\, 0,\, h_z]^\top$ [7, Sec. 3.4, Eqs. (45)--(46)]. The
magnetic term is a full cross product, so it also corrects roll and pitch. A
sensor reading of exactly zero contributes no term.

Complementary filters are less optimal than the MEKF but require no covariance
propagation, no matrix inversions, and have deterministic constant-time
execution &mdash; making them suitable for resource-constrained embedded systems.

## References

[1] T. D. Barfoot, "State Estimation for Robotics," 2nd ed., Cambridge
University Press, 2024.

[2] K. M. Lynch and F. C. Park, "Modern Robotics: Mechanics, Planning, and
Control," Cambridge University Press, 2017.

[3] F. L. Markley, "Attitude Error Representations for Kalman Filtering,"
Journal of Guidance, Control, and Dynamics, vol. 26, no. 2, pp. 311--317,
2003.

[4] S. Hauberg, F. Lauze, and K. S. Pedersen, "Unscented Kalman Filtering on
Riemannian Manifolds," Journal of Mathematical Imaging and Vision, vol.
46, no. 1, pp. 103--120, 2013.

[5] R. Mahony, T. Hamel, and J.-M. Pflimlin, "Nonlinear Complementary
Filters on the Special Orthogonal Group," IEEE Transactions on Automatic
Control, vol. 53, no. 5, pp. 1203--1218, 2008.

[6] J. Sola, J. Deray, and D. Atchuthan, "A Micro Lie Theory for State
Estimation in Robotics," arXiv:1812.01537v9, 2021.

[7] S. O. H. Madgwick, "An Efficient Orientation Filter for Inertial and
Inertial/Magnetic Sensor Arrays," report, 2010.
