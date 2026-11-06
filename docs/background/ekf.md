# Extended Kalman Filter

The Extended Kalman Filter (EKF) extends the linear Kalman filter to nonlinear
systems by linearizing the dynamics and measurement models around the current
state estimate at each time step. It is the most widely used nonlinear state
estimator due to its simplicity and efficiency for mildly nonlinear systems
[1, Ch. 13, p. 396].

The EKF was developed in the 1960s for aerospace navigation and remains the
workhorse estimator in many applications including GPS, inertial navigation,
and robotic localization.

## Nonlinear System Model

The EKF addresses discrete nonlinear systems of the form
[1, Sec. 13.2.3, Eq. (13.38), p. 407]:

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
evaluated at each time step [1, Sec. 13.2.3, Eqs. (13.46) and (13.48), p. 409]:

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
prediction uses the linearized dynamics [1, Sec. 13.2.3, Eq. (13.47), p. 409]:

$$
\hat{x}_{k+1|k} = f(\hat{x}_{k|k}, u_k)
$$

$$
P_{k+1|k} = F_k \, P_{k|k} \, F_k^\top + Q
$$

Note that the state propagation applies $f$ directly (not its linearization),
preserving accuracy in the mean. Only the covariance propagation uses the
Jacobian approximation.

## EKF Update

The measurement update follows the standard Kalman structure with the
nonlinear measurement function and its Jacobian
[1, Sec. 13.2.3, Eq. (13.49), p. 409]:

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
the predicted state, not a linearized version.

ctrlpp's `ekf` forms the covariance update in the Joseph form,
$(I - K_k H_k) P_{k|k-1} (I - K_k H_k)^\top + K_k R K_k^\top$, which equals
the short form above for this $K_k$ in exact arithmetic and stays positive
semidefinite for any gain, and it symmetrizes both the predicted and the
updated covariance.

## Jacobian Computation

The EKF requires Jacobians at every time step. Two approaches exist:

### Analytical Jacobians

The user provides closed-form expressions for $F_k$ and $H_k$. This is the
most efficient approach when the derivatives are tractable. For example,
if $f(x) = [x_1 x_2, \; x_1 + x_2^2]^\top$, then

$$
F = \begin{bmatrix} x_2 & x_1 \\ 1 & 2 x_2 \end{bmatrix}
$$

### Numerical Jacobians

When analytical derivatives are unavailable or impractical, finite-difference
approximations compute the Jacobians:

$$
F_{ij} \approx \frac{f_i(x + \delta e_j) - f_i(x - \delta e_j)}{2\delta}
$$

where $e_j$ is the $j$-th standard basis vector and $\delta$ is a small
perturbation. Central differences provide second-order accuracy in $\delta$.
ctrlpp perturbs component $j$ by
$\delta_j = \epsilon^{1/3} \max(1, |x_j|)$, rounded to the representable
difference $(x_j + \delta_j) - x_j$; the cube root of machine epsilon
balances the $O(\delta^2)$ truncation error against the $O(\epsilon / \delta)$
rounding error. The factor $\epsilon^{1/3}$ is `ekf_config::numerical_eps`.

## Limitations

The EKF has several well-known limitations
[1, Sec. 13.5, p. 425; Ch. 14, p. 433]:

1. **First-order accuracy**: the linearization captures only the first
   derivative, introducing errors proportional to the neglected higher-order
   terms. For strongly nonlinear systems, these errors can cause divergence.

2. **Jacobian requirement**: computing Jacobians can be error-prone and
   computationally expensive for high-dimensional systems.

3. **Non-guaranteed convergence**: unlike the linear Kalman filter, the EKF
   has no general convergence guarantees. Initialization errors, model
   mismatch, or strong nonlinearities can cause the filter to diverge.

4. **Gaussian assumption**: the EKF assumes the posterior distribution
   remains approximately Gaussian. Multi-modal or heavily skewed
   distributions are poorly represented.

## Iterated EKF

The iterated EKF (IEKF) re-linearizes the measurement model around the
updated state estimate and repeats the update step until convergence
[1, Sec. 13.3.1, pp. 410--413]. This is equivalent to performing a
Gauss-Newton optimization on the measurement update and can improve
accuracy for strongly nonlinear measurement functions. ctrlpp's `ekf` does not
iterate: it linearizes the measurement once per update.

## Consistency Monitoring

The innovation sequence $\tilde{y}_k$ provides a diagnostic for filter
health. For a well-tuned EKF:

- The normalized innovation squared (NIS):
  $\epsilon_k = \tilde{y}_k^\top S_k^{-1} \tilde{y}_k$ should follow a
  $\chi^2$ distribution with $m$ degrees of freedom
- The normalized estimation error squared (NEES):
  $\epsilon_k = \tilde{x}_k^\top P_{k|k}^{-1} \tilde{x}_k$ should follow
  $\chi^2$ with $n$ degrees of freedom (when true state is available)

`ekf::nis()` returns the NIS of the last update; the NEES needs the true state
and is left to the caller.

Persistent departures from these distributions indicate model mismatch or
incorrect noise tuning.

## References

[1] D. Simon, "Optimal State Estimation: Kalman, H-Infinity, and Nonlinear
Approaches," Wiley, 2006.
