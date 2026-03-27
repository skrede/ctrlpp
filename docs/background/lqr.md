# Linear Quadratic Regulator

The Linear Quadratic Regulator (LQR) is the optimal state-feedback controller
for linear systems with a quadratic cost function. It provides a systematic
method for computing gain matrices that balance state regulation against
control effort, with guaranteed stability margins and a closed-form solution
via the Riccati equation [1, Ch. 12, pp. 671--720].

LQR is the foundation of modern optimal control and forms the basis for many
advanced control strategies including LQG (LQR + Kalman filter), MPC terminal
cost design, and integral-action state feedback.

## Optimal Control Problem

Consider a discrete linear time-invariant system [1, Sec. 12.8, pp. 709--714]:

$$
x_{k+1} = A \, x_k + B \, u_k
$$

The LQR problem seeks the control sequence $\{u_0, u_1, \ldots\}$ that
minimises the infinite-horizon quadratic cost [2, Sec. 3.2, pp. 43--50]:

$$
J = \sum_{k=0}^{\infty} \left( x_k^\top Q \, x_k + u_k^\top R \, u_k \right)
$$

where:

- $Q \succeq 0$ (positive semi-definite) is the state weight matrix, penalising
  deviation from the origin
- $R \succ 0$ (positive definite) is the input weight matrix, penalising
  control effort

The weight matrices $Q$ and $R$ encode the designer's trade-off: larger $Q$
relative to $R$ produces aggressive control with tight regulation, while
larger $R$ relative to $Q$ produces conservative control that limits actuator
usage.

## Discrete Algebraic Riccati Equation

The optimal cost-to-go matrix $P$ satisfies the Discrete Algebraic Riccati
Equation (DARE) [1, Sec. 12.8, p. 711]:

$$
P = A^\top P \, A - A^\top P \, B \, (R + B^\top P \, B)^{-1} B^\top P \, A + Q
$$

This is a matrix fixed-point equation that can be solved iteratively or via
eigendecomposition methods. The solution exists and is unique when the system
$(A, B)$ is stabilisable and $(A, Q^{1/2})$ is detectable
[2, Sec. 3.3, pp. 51--55].

### Schur Method

The DARE can be solved reliably using the generalised Schur (QZ) decomposition
of the symplectic pencil [3, pp. 913--921]:

$$
\begin{bmatrix} A + B R^{-1} B^\top (A^\top)^{-1} Q & -B R^{-1} B^\top (A^\top)^{-1} \\ -(A^\top)^{-1} Q & (A^\top)^{-1} \end{bmatrix}
$$

The stable eigenspace of this pencil yields the solution $P$. This approach
avoids the numerical issues of iterative methods and handles ill-conditioned
problems robustly.

## Optimal Gain

Given the DARE solution $P$, the optimal state-feedback gain is
[1, Sec. 12.8, p. 712]:

$$
K = (R + B^\top P \, B)^{-1} B^\top P \, A
$$

The optimal control law is then:

$$
u_k = -K \, x_k
$$

This is a static linear state-feedback law: the gain matrix $K$ is computed
offline and applied at each time step with no online optimisation required.

## Stability Properties

The LQR has strong stability guarantees [2, Sec. 3.4, pp. 56--60]:

- **Closed-loop stability**: all eigenvalues of $(A - BK)$ lie strictly
  inside the unit circle (for discrete-time systems).
- **Gain margin**: for continuous-time LQR, the gain margin is at least
  $[1/2, \infty)$ and the phase margin is at least $60^\circ$.
- **Robustness**: the LQR provides inherent robustness to model uncertainty,
  though the discrete-time margins are generally less generous than
  continuous-time.

## Finite-Horizon LQR

For a finite horizon of $N$ steps, the cost function is
[1, Sec. 12.8, pp. 710--711]:

$$
J = \sum_{k=0}^{N-1} \left( x_k^\top Q \, x_k + u_k^\top R \, u_k \right) + x_N^\top Q_f \, x_N
$$

where $Q_f$ is the terminal cost weight. The optimal gain varies with time
and is computed by iterating the discrete Riccati recursion backward from
$P_N = Q_f$ [2, Sec. 3.1, pp. 40--43]:

$$
P_k = A^\top P_{k+1} A - A^\top P_{k+1} B (R + B^\top P_{k+1} B)^{-1} B^\top P_{k+1} A + Q
$$

$$
K_k = (R + B^\top P_{k+1} B)^{-1} B^\top P_{k+1} A
$$

As $N \to \infty$, the time-varying gains converge to the steady-state
infinite-horizon gain.

## LQR with Integral Action

Standard LQR regulates states to the origin. To track non-zero references
and reject constant disturbances, the state vector is augmented with integral
states [1, Sec. 12.9, pp. 715--718]:

$$
x_{\text{aug}} = \begin{bmatrix} x_k \\ \xi_k \end{bmatrix}, \qquad
\xi_{k+1} = \xi_k + T_s (r_k - C x_k)
$$

The augmented system is:

$$
x_{\text{aug}, k+1} = \begin{bmatrix} A & 0 \\ -T_s C & I \end{bmatrix} x_{\text{aug}, k} + \begin{bmatrix} B \\ 0 \end{bmatrix} u_k + \begin{bmatrix} 0 \\ T_s I \end{bmatrix} r_k
$$

Solving the DARE for this augmented system yields gains for both the state
feedback and integral action components. The integral states provide the same
benefit as the I-term in PID: zero steady-state error for step references and
step disturbances.

## Weight Selection

The choice of $Q$ and $R$ determines the controller characteristics
[1, Sec. 12.8, pp. 713--714]:

- **Bryson's rule**: set $Q_{ii} = 1 / x_{i,\max}^2$ and
  $R_{jj} = 1 / u_{j,\max}^2$, where $x_{i,\max}$ and $u_{j,\max}$ are the
  maximum acceptable values. This normalises the cost terms.
- **Iterative tuning**: start with identity matrices and adjust based on
  closed-loop simulation.
- **Cross-weight $N$**: the generalised cost
  $J = \sum (x^\top Q x + u^\top R u + 2 x^\top N u)$ allows penalising
  state-input correlations.

## References

[1] N. S. Nise, "Control Systems Engineering," 7th ed., Wiley, 2015.

[2] D. Simon, "Optimal State Estimation: Kalman, H-Infinity, and Nonlinear
Approaches," Wiley, 2006.

[3] A. J. Laub, "A Schur Method for Solving Algebraic Riccati Equations,"
IEEE Transactions on Automatic Control, vol. 24, no. 6, pp. 913--921, 1979.

[4] B. D. O. Anderson and J. B. Moore, "Optimal Control: Linear Quadratic
Methods," Dover, 1990.
