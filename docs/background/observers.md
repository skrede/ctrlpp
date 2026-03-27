# State Observers

A state observer reconstructs the internal state of a dynamic system from
its measured inputs and outputs. When only a subset of the state variables
can be measured directly, observers provide estimates of the full state vector
for use in state-feedback control laws such as LQR
[1, Sec. 12.5, pp. 685--700].

The observer concept was introduced by Luenberger in 1964 and formalised in
his 1971 paper [2], establishing a deterministic dual to the Kalman filter.

## The Observation Problem

Given a discrete linear time-invariant system
[1, Sec. 12.5, pp. 685--688]:

$$
x_{k+1} = A \, x_k + B \, u_k
$$

$$
y_k = C \, x_k
$$

where $x_k \in \mathbb{R}^n$ is the state, $u_k \in \mathbb{R}^p$ is the
known input, and $y_k \in \mathbb{R}^m$ is the measured output. The goal is
to construct an estimate $\hat{x}_k$ that converges to the true state
$x_k$ as $k \to \infty$, using only the known sequences $\{u_k\}$ and
$\{y_k\}$.

## Luenberger Observer

The Luenberger observer augments a copy of the system model with a correction
term driven by the output estimation error
[1, Sec. 12.5, pp. 688--692]:

$$
\hat{x}_{k+1} = A \, \hat{x}_k + B \, u_k + L \, (y_k - C \, \hat{x}_k)
$$

where $L \in \mathbb{R}^{n \times m}$ is the observer gain matrix. The
estimation error $e_k = x_k - \hat{x}_k$ evolves according to:

$$
e_{k+1} = (A - L C) \, e_k
$$

The error converges to zero if and only if all eigenvalues of $(A - LC)$ lie
strictly inside the unit circle. The observer gain $L$ determines how quickly
the estimate converges.

## Observability

The observer can reconstruct the full state only if the system is observable.
The observability matrix [1, Sec. 12.3, pp. 683--685]:

$$
\mathcal{O} = \begin{bmatrix} C \\ CA \\ CA^2 \\ \vdots \\ CA^{n-1} \end{bmatrix}
$$

must have rank $n$ (full column rank). If $\text{rank}(\mathcal{O}) < n$,
only the observable subspace can be estimated. The unobservable modes evolve
according to the open-loop dynamics and cannot be influenced by the observer.

### Detectability

A weaker condition than full observability is detectability: all unobservable
modes must be stable (eigenvalues inside the unit circle). This is sufficient
for the observer error to converge, since the unobservable modes decay
naturally [1, Sec. 12.5, p. 690].

## Observer Pole Placement

The observer gain $L$ is chosen to place the eigenvalues of $(A - LC)$ at
desired locations. This is the dual of the state-feedback pole placement
problem [1, Sec. 12.5, pp. 692--696].

### Ackermann's Formula

For single-output systems ($m = 1$), Ackermann's formula provides a direct
computation [1, Sec. 12.5, pp. 693--694]:

$$
L = \alpha_d(A) \, \mathcal{O}^{-1} \begin{bmatrix} 0 \\ 0 \\ \vdots \\ 1 \end{bmatrix}
$$

where $\alpha_d(s) = \prod_{i=1}^{n} (s - p_i)$ is the desired characteristic
polynomial evaluated at $A$.

### Multi-Output Placement

For multi-output systems, more general algorithms such as the method of
Kautsky, Nichols, and Van Dooren [3, pp. 1129--1155] are used. This method
iteratively selects eigenvectors to minimise sensitivity while placing
eigenvalues at the specified locations.

## Observer Design Guidelines

The choice of observer poles involves a trade-off
[1, Sec. 12.5, pp. 696--698]:

- **Fast poles** (close to origin): rapid convergence but high sensitivity to
  measurement noise and modelling errors.
- **Slow poles** (close to unit circle): smooth estimates but slow tracking
  of actual state changes.

A common heuristic is to place observer poles 2--6 times faster than the
closed-loop controller poles [1, Sec. 12.5, p. 697]. This ensures the
observer converges before the controller needs accurate state information.

## Separation Principle

The separation principle states that the state-feedback gain $K$ and the
observer gain $L$ can be designed independently: the closed-loop eigenvalues
of the combined controller-observer system are the union of the controller
eigenvalues and the observer eigenvalues
[1, Sec. 12.6, pp. 699--704].

This means:

1. Design $K$ assuming perfect state measurement (state-feedback design)
2. Design $L$ independently (observer design)
3. The combined system $u_k = -K \hat{x}_k$ with observer has eigenvalues
   $\text{eig}(A - BK) \cup \text{eig}(A - LC)$

The separation principle holds for linear systems but does not extend to
nonlinear systems in general.

## Relationship to Kalman Filter

The Kalman filter is a stochastic observer that computes the gain $L = K_k$
(the Kalman gain) optimally from noise statistics $(Q, R)$. The Luenberger
observer is its deterministic counterpart where $L$ is chosen by pole
placement rather than optimisation [1, Sec. 12.7, pp. 705--708].

When the noise covariances are known, the Kalman filter is preferred. When
they are unknown or the designer wants direct control over convergence speed,
the Luenberger observer with pole placement is more transparent.

## Reduced-Order Observers

A full-order observer estimates all $n$ states, including those that are
directly measured. A reduced-order (Luenberger) observer estimates only the
$n - m$ unmeasured states, using the measured outputs directly. This reduces
computation but increases design complexity
[1, Sec. 12.5, pp. 698--699].

## References

[1] N. S. Nise, "Control Systems Engineering," 7th ed., Wiley, 2015.

[2] D. G. Luenberger, "An Introduction to Observers," IEEE Transactions on
Automatic Control, vol. 16, no. 6, pp. 596--602, 1971.

[3] J. Kautsky, N. K. Nichols, and P. Van Dooren, "Robust Pole Assignment in
Linear State Feedback," International Journal of Control, vol. 41, no. 5,
pp. 1129--1155, 1985.

[4] D. Simon, "Optimal State Estimation: Kalman, H-Infinity, and Nonlinear
Approaches," Wiley, 2006.
