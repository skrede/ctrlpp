# SO(3) and Lie Groups

The special orthogonal group SO(3) is the group of all rotations in
three-dimensional space. It forms a Lie group: a smooth manifold that is also
a group, meaning rotations can be composed and inverted smoothly. Understanding
SO(3) is essential for attitude estimation, robotic manipulation, and any
application involving 3D orientation [1, Ch. 7, pp. 211--266].

This page covers rotation representations, the Lie group structure, the
exponential and logarithmic maps, and the Adjoint representation.

## Rotation Matrices

A rotation matrix $R \in \mathbb{R}^{3 \times 3}$ satisfies the orthogonality
constraints [1, Sec. 7.1.1, pp. 213--215]:

$$
R^\top R = I, \qquad \det(R) = +1
$$

The set of all such matrices forms SO(3):

$$
\text{SO}(3) = \{ R \in \mathbb{R}^{3 \times 3} : R^\top R = I, \; \det(R) = 1 \}
$$

SO(3) is a 3-dimensional manifold embedded in $\mathbb{R}^{3 \times 3}$. It
has 9 entries but only 3 degrees of freedom, constrained by 6 orthogonality
conditions.

### Rotation Matrix Composition

Rotations compose by matrix multiplication. If $R_1$ rotates frame $A$ to $B$
and $R_2$ rotates frame $B$ to $C$, then $R_2 R_1$ rotates $A$ to $C$
[1, Sec. 7.1.1, p. 214].

The inverse of a rotation is its transpose: $R^{-1} = R^\top$.

## Unit Quaternions

A unit quaternion $q = [q_w, q_x, q_y, q_z]^\top$ with $\lVert q \rVert = 1$
provides a compact, singularity-free representation of rotations. Quaternions
form a double cover of SO(3): both $q$ and $-q$ represent the same rotation
[1, Sec. 7.1.3, pp. 218--222].

### Hamilton Product

The Hamilton product composes two rotations [1, eq. (7.12), p. 219]:

$$
p \otimes q = \begin{bmatrix}
p_w q_w - p_x q_x - p_y q_y - p_z q_z \\
p_w q_x + p_x q_w + p_y q_z - p_z q_y \\
p_w q_y - p_x q_z + p_y q_w + p_z q_x \\
p_w q_z + p_x q_y - p_y q_x + p_z q_w
\end{bmatrix}
$$

The inverse quaternion is the conjugate: $q^{-1} = [q_w, -q_x, -q_y, -q_z]^\top$.

### Quaternion to Rotation Matrix

The equivalent rotation matrix is [1, eq. (7.14), p. 220]:

$$
R(q) = (q_w^2 - \lVert q_{xyz} \rVert^2) I + 2 q_{xyz} q_{xyz}^\top + 2 q_w [q_{xyz}]_\times
$$

## Axis-Angle Representation

A rotation can be described by an axis $\hat{e} \in \mathbb{R}^3$ with
$\lVert \hat{e} \rVert = 1$ and an angle $\theta \in [0, \pi]$. The rotation
vector $\phi = \theta \hat{e} \in \mathbb{R}^3$ is a minimal (3-parameter)
representation [1, Sec. 7.1.2, pp. 215--218].

The rotation vector is undefined when $\theta = 0$ (the axis direction is
arbitrary for the identity rotation). This singularity at the identity is
unavoidable for any 3-parameter representation (a topological constraint).

## The Lie Algebra so(3)

The Lie algebra of SO(3), denoted $\mathfrak{so}(3)$, is the tangent space at
the identity element. It consists of all $3 \times 3$ skew-symmetric matrices
[1, Sec. 7.2.1, pp. 223--225]:

$$
\mathfrak{so}(3) = \{ [\phi]_\times : \phi \in \mathbb{R}^3 \}
$$

where the hat (wedge) operator maps a 3-vector to its skew-symmetric matrix:

$$
[\phi]_\times = \begin{bmatrix} 0 & -\phi_3 & \phi_2 \\ \phi_3 & 0 & -\phi_1 \\ -\phi_2 & \phi_1 & 0 \end{bmatrix}
$$

The inverse (vee) operator recovers $\phi$ from $[\phi]_\times$.

## Exponential Map

The exponential map takes an element of the Lie algebra (a tangent vector)
to the Lie group (a rotation). For $\phi \in \mathbb{R}^3$ with
$\theta = \lVert \phi \rVert$ [1, eq. (7.18), p. 229]:

$$
\exp([\phi]_\times) = I + \frac{\sin \theta}{\theta} [\phi]_\times + \frac{1 - \cos \theta}{\theta^2} [\phi]_\times^2
$$

This is Rodrigues' rotation formula. Near $\theta = 0$, the sinc-type
coefficients are evaluated using Taylor series to avoid division by zero:

$$
\frac{\sin \theta}{\theta} \approx 1 - \frac{\theta^2}{6} + \frac{\theta^4}{120}, \qquad
\frac{1 - \cos \theta}{\theta^2} \approx \frac{1}{2} - \frac{\theta^2}{24} + \frac{\theta^4}{720}
$$

For quaternions, the exponential map is [1, eq. (7.20), p. 230]:

$$
\exp(\phi) = \begin{bmatrix} \cos(\theta/2) \\ \frac{\sin(\theta/2)}{\theta} \phi \end{bmatrix}
$$

## Logarithmic Map

The logarithmic map is the inverse of the exponential map, taking a rotation
back to the Lie algebra [1, eq. (7.22), p. 231]:

For a rotation matrix $R$ with rotation angle
$\theta = \arccos\!\left(\frac{\text{tr}(R) - 1}{2}\right)$:

$$
\log(R) = \frac{\theta}{2 \sin \theta} (R - R^\top)
$$

For a quaternion $q$:

$$
\log(q) = \frac{2 \arccos(q_w)}{\lVert q_{xyz} \rVert} \, q_{xyz}
$$

with sinc-expansion near $\theta = 0$ for numerical stability.

## Adjoint Representation

The Adjoint representation of SO(3) describes how the Lie algebra transforms
under conjugation by group elements
[1, Sec. 7.3, pp. 234--238]:

$$
\text{Ad}_R : \mathfrak{so}(3) \to \mathfrak{so}(3), \qquad \text{Ad}_R(\phi) = R \phi
$$

For SO(3), the Adjoint representation is simply the rotation matrix itself
applied to tangent vectors. This is used extensively in:

- Transforming angular velocities between frames
- Propagating uncertainty on the manifold
- Computing Jacobians of composed rotations

## Left and Right Jacobians

The left Jacobian of SO(3) relates perturbations in the Lie algebra to
perturbations on the group [1, Sec. 7.4, pp. 240--245]:

$$
J_l(\phi) = I + \frac{1 - \cos \theta}{\theta^2} [\phi]_\times + \frac{\theta - \sin \theta}{\theta^3} [\phi]_\times^2
$$

The right Jacobian is $J_r(\phi) = J_l(-\phi) = J_l(\phi)^\top$. These
Jacobians appear in:

- Error-state Kalman filter (MEKF) covariance propagation
- Optimisation on manifolds (pose graph SLAM)
- Interpolation between rotations

## Composition and Perturbation

For estimation and control on SO(3), it is common to represent uncertainty as
a small perturbation $\delta\phi \in \mathbb{R}^3$ in the tangent space
[1, Sec. 7.5, pp. 245--250]:

$$
R_{\text{true}} = R_{\text{nominal}} \cdot \exp([\delta\phi]_\times)
$$

This separates the manifold structure (handled by the nominal rotation) from
the linear perturbation (handled by standard linear algebra). The MEKF and
manifold UKF both use this decomposition.

## References

[1] T. D. Barfoot, "State Estimation for Robotics," 2nd ed., Cambridge
University Press, 2025.

[2] K. M. Lynch and F. C. Park, "Modern Robotics: Mechanics, Planning, and
Control," Cambridge University Press, 2017.

[3] J. Sola, J. Deray, and D. Atchuthan, "A Micro Lie Theory for State
Estimation in Robotics," arXiv:1812.01537, 2018.
