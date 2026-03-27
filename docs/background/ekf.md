# Extended Kalman Filter

The Extended Kalman Filter (EKF) extends the linear Kalman filter to nonlinear
systems by linearising the dynamics and measurement models around the current
state estimate at each time step. It is the most widely used nonlinear state
estimator due to its simplicity and efficiency for mildly nonlinear systems
[1, Ch. 13, pp. 393--420].

The EKF was developed in the 1960s for aerospace navigation and remains the
workhorse estimator in many applications including GPS, inertial navigation,
and robotic localisation.

## Nonlinear System Model

The EKF addresses discrete nonlinear systems of the form
[1, Sec. 13.1, pp. 394--396]:

$$
x_{k+1} = f(x_k, u_k) + w_k, \qquad w_k \sim \mathcal{N}(0, Q)
$$

$$
z_k = h(x_k) + v_k, \qquad v_k \sim \mathcal{N}(0, R)
$$

where $f$ is the nonlinear process model, $h$ is the nonlinear measurement
model, and the noise terms remain additive and Gaussian. Unlike the linear
Kalman filter, the standard Kalman equations cannot be applied directly because
the transformations of Gaussian random variables through nonlinear functions
are no longer Gaussian in general.

## First-Order Taylor Expansion

The EKF approximates the nonlinear functions by their first-order Taylor
series expansions around the current estimate. The Jacobian matrices are
evaluated at each time step [1, Sec. 13.1, pp. 397--399]:

$$
F_k = \left. \frac{\partial f}{\partial x} \right|_{\hat{x}_{k|k}, u_k} \in \mathbb{R}^{n \times n}
$$

$$
H_k = \left. \frac{\partial h}{\partial x} \right|_{\hat{x}_{k|k-1}} \in \mathbb{R}^{m \times n}
$$

These Jacobians replace the constant system matrices $A$ and $C$ in the
standard Kalman filter equations. The quality of the EKF depends directly on
how well the first-order approximation captures the true nonlinearity.

## EKF Prediction

The state prediction uses the full nonlinear model, while the covariance
prediction uses the linearised dynamics [1, Sec. 13.2, pp. 400--403]:

$$
\hat{x}_{k+1|k} = f(\hat{x}_{k|k}, u_k)
$$

$$
P_{k+1|k} = F_k \, P_{k|k} \, F_k^\top + Q
$$

Note that the state propagation applies $f$ directly (not its linearisation),
preserving accuracy in the mean. Only the covariance propagation uses the
Jacobian approximation.

## EKF Update

The measurement update follows the standard Kalman structure with the
nonlinear measurement function and its Jacobian
[1, Sec. 13.2, pp. 403--405]:

$$
\tilde{y}_k = z_k - h(\hat{x}_{k|k-1})
$$

$$
S_k = H_k \, P_{k|k-1} \, H_k^\top + R
$$

$$
K_k = P_{k|k-1} \, H_k^\top \, S_k^{-1}
$$

$$
\hat{x}_{k|k} = \hat{x}_{k|k-1} + K_k \, \tilde{y}_k
$$

$$
P_{k|k} = (I - K_k \, H_k) \, P_{k|k-1}
$$

The innovation $\tilde{y}_k$ uses the nonlinear function $h$ evaluated at
the predicted state, not a linearised version.

## Jacobian Computation

The EKF requires Jacobians at every time step. Two approaches exist:

### Analytical Jacobians

The user provides closed-form expressions for $F_k$ and $H_k$. This is the
most efficient approach when the derivatives are tractable. For example,
if $f(x) = [x_1 x_2, \; x_1 + x_2^2]^\top$, then
[1, Sec. 13.3, pp. 406--408]:

$$
F = \begin{bmatrix} x_2 & x_1 \\ 1 & 2 x_2 \end{bmatrix}
$$

### Numerical Jacobians

When analytical derivatives are unavailable or impractical, finite-difference
approximations compute the Jacobians [1, Sec. 13.3, pp. 408--410]:

$$
F_{ij} \approx \frac{f_i(x + \delta e_j) - f_i(x - \delta e_j)}{2\delta}
$$

where $e_j$ is the $j$-th standard basis vector and $\delta$ is a small
perturbation. Central differences provide second-order accuracy in $\delta$.

## Limitations

The EKF has several well-known limitations
[1, Sec. 13.4, pp. 411--415]:

1. **First-order accuracy**: the linearisation captures only the first
   derivative, introducing errors proportional to the neglected higher-order
   terms. For strongly nonlinear systems, these errors can cause divergence.

2. **Jacobian requirement**: computing Jacobians can be error-prone and
   computationally expensive for high-dimensional systems.

3. **Non-guaranteed convergence**: unlike the linear Kalman filter, the EKF
   has no general convergence guarantees. Initialisation errors, model
   mismatch, or strong nonlinearities can cause the filter to diverge.

4. **Gaussian assumption**: the EKF assumes the posterior distribution
   remains approximately Gaussian. Multi-modal or heavily skewed
   distributions are poorly represented.

## Iterated EKF

The iterated EKF (IEKF) re-linearises the measurement model around the
updated state estimate and repeats the update step until convergence
[1, Sec. 13.5, pp. 416--418]. This is equivalent to performing a
Gauss-Newton optimisation on the measurement update and can improve
accuracy for strongly nonlinear measurement functions.

## Consistency Monitoring

The innovation sequence $\tilde{y}_k$ provides a diagnostic for filter
health. For a well-tuned EKF [2, Sec. 8.4, pp. 232--235]:

- The normalised innovation squared (NIS):
  $\epsilon_k = \tilde{y}_k^\top S_k^{-1} \tilde{y}_k$ should follow a
  $\chi^2$ distribution with $m$ degrees of freedom
- The normalised estimation error squared (NEES):
  $\epsilon_k = \tilde{x}_k^\top P_{k|k}^{-1} \tilde{x}_k$ should follow
  $\chi^2$ with $n$ degrees of freedom (when true state is available)

Persistent departures from these distributions indicate model mismatch or
incorrect noise tuning.

## References

[1] D. Simon, "Optimal State Estimation: Kalman, H-Infinity, and Nonlinear
Approaches," Wiley, 2006.

[2] S. Sarkka and L. Svensson, "Bayesian Filtering and Smoothing," 2nd ed.,
Cambridge University Press, 2023.

[3] S. Thrun, W. Burgard, and D. Fox, "Probabilistic Robotics," MIT Press,
2006.
