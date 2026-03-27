# Moving Horizon Estimation

Moving Horizon Estimation (MHE) is the estimation dual of Model Predictive
Control. Instead of optimising future inputs over a prediction horizon, MHE
optimises past state estimates over a fixed window of recent measurements. It
provides optimal state estimation for constrained systems where the Kalman
filter cannot enforce physical bounds on states
[1, Sec. 4, pp. 248--252].

MHE was formalised by Rao, Rawlings, and Mayne [1] and shares the same
receding-horizon philosophy as MPC: solve a finite-window optimisation at each
step, use only the most recent estimate, and slide the window forward.

## The MHE Cost Function

At time $k$, MHE estimates the state trajectory over the most recent $N$ time
steps. For a linear system $x_{i+1} = A x_i + B u_i + w_i$ with measurements
$z_i = C x_i + v_i$, the decision variables are the initial state $x_{k-N}$
and the process noise sequence $\{w_{k-N}, \ldots, w_{k-1}\}$
[1, Sec. 3, pp. 249--253]:

$$
\min_{x_{k-N}, \, w_{k-N:k-1}} \;
\lVert x_{k-N} - \bar{x}_{k-N} \rVert^2_{P_{\text{arr}}^{-1}}
+ \sum_{i=k-N}^{k-1} \lVert w_i \rVert^2_{Q^{-1}}
+ \sum_{i=k-N}^{k} \lVert z_i - C x_i \rVert^2_{R^{-1}}
$$

$$
\text{subject to} \quad x_{i+1} = A x_i + B u_i + w_i
$$

$$
x_{\min} \le x_i \le x_{\max}
$$

where the weighted norm is $\lVert a \rVert^2_M = a^\top M \, a$.

The three terms balance:
1. **Arrival cost**: consistency with prior information before the window
2. **Process noise**: penalising deviations from the model
3. **Measurement fit**: consistency with observations

## Arrival Cost

The first term $\lVert x_{k-N} - \bar{x}_{k-N} \rVert^2_{P_{\text{arr}}^{-1}}$
is the arrival cost. It summarises all information from measurements before
the estimation window into a single quadratic penalty
[1, Sec. 3.1, pp. 250--251].

As the window slides forward, the arrival cost must be updated:

- **Prior mean** $\bar{x}_{k-N}$: the best estimate of the state at the
  window boundary, typically from a companion filter
- **Arrival cost weighting** $P_{\text{arr}}$: encodes the uncertainty of
  that estimate

### Exact Update (Linear Systems)

For linear systems without constraints, the exact arrival cost equals the
Kalman filter covariance: $\bar{x}_{k-N} = \hat{x}_{k-N}^{\text{KF}}$ and
$P_{\text{arr}} = P_{k-N}^{\text{KF}}$ [1, Sec. 3.2, pp. 251--253].

### Approximate Update (Nonlinear Systems)

For nonlinear systems, an EKF approximation provides the arrival cost. A
companion EKF runs alongside the MHE, and its state estimate and covariance
are used as $\bar{x}$ and $P_{\text{arr}}$ at each window shift
[2, Sec. 4.3, pp. 403--405].

## Relationship to the Kalman Filter

For linear systems without state constraints, MHE with infinite window length
recovers the Kalman filter solution exactly. With finite window, MHE
approximates the Kalman filter but can additionally enforce
[1, Sec. 4.1, pp. 253--255]:

$$
x_{\min} \le x_i \le x_{\max}
$$

This is the key advantage of MHE: state constraints are handled naturally
within the optimisation. Physical constraints such as non-negative
concentrations, bounded temperatures, or positive definite covariances
can be directly enforced.

## QP Formulation

For linear dynamics, the MHE problem is a quadratic program. The structure
is similar to linear MPC with the time direction reversed
[3, Sec. 4.2, pp. 400--403]:

$$
\min_z \; \frac{1}{2} z^\top H z + q^\top z
$$

$$
\text{subject to} \quad A_{\text{eq}} z = b_{\text{eq}}
$$

$$
l \le z \le u
$$

where $z$ collects the state trajectory and noise variables. The QP has the
same block-banded sparsity pattern as linear MPC, enabling efficient solution
with the same solvers.

## Convergence and Stability

Rao, Rawlings, and Mayne [1, Sec. 5, pp. 255--258] proved that MHE
converges under the following conditions:

1. The system is detectable (all unobservable modes are stable)
2. The arrival cost is updated consistently (EKF approximation or exact)
3. The estimation window $N$ is sufficiently long

The estimation error is bounded and converges to zero as the window length
increases. In practice, window lengths of 10--50 steps are typical, depending
on the system dynamics and measurement rate.

## Tuning Guidelines

- **Window length $N$**: longer windows use more data but increase
  computation. Start with $N = 10$--$20$.
- **$Q$ and $R$**: same interpretation as in the Kalman filter. $Q$
  reflects model uncertainty, $R$ reflects measurement noise.
- **State bounds**: should reflect true physical constraints, not
  artificially tight bounds (which slow convergence).

## References

[1] C. V. Rao, J. B. Rawlings, and D. Q. Mayne, "Constrained State Estimation
for Nonlinear Discrete-Time Systems: Stability and Moving Horizon
Approximations," IEEE Transactions on Automatic Control, vol. 48, no. 2,
pp. 246--258, 2003.

[2] M. Diehl, H. J. Ferreau, and N. Haverbeke, "Efficient Numerical Methods
for Nonlinear MPC and Moving Horizon Estimation," in Nonlinear Model
Predictive Control, LNCIS vol. 384, Springer, 2009, pp. 391--417.

[3] J. B. Rawlings, D. Q. Mayne, and M. Diehl, "Model Predictive Control:
Theory, Computation, and Design," 2nd ed., Nob Hill Publishing, 2017.

[4] D. Simon, "Optimal State Estimation: Kalman, H-Infinity, and Nonlinear
Approaches," Wiley, 2006.
