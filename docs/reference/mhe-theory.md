# MHE Theory

Moving Horizon Estimation (MHE) is the estimation dual of Model Predictive
Control. Instead of optimising future inputs over a prediction horizon, MHE
optimises past state estimates over a fixed window of recent measurements.
It handles state constraints and non-Gaussian noise naturally.

## Key Concepts

### Finite-Window Estimation

MHE estimates the state trajectory over the most recent N time steps by
solving an optimisation problem:

```
minimise  ||x[k-N] - x_bar||^2_{P_bar^{-1}}
        + sum_{i=k-N}^{k-1} ||w[i]||^2_{Q^{-1}}
        + sum_{i=k-N}^{k} ||v[i]||^2_{R^{-1}}

subject to  x[i+1] = f(x[i], u[i]) + w[i]
            z[i]   = h(x[i]) + v[i]
            x[i] in X   (state constraints)
```

The decision variables are the initial state x[k-N] and the process noise
sequence w[k-N:k-1].

### Arrival Cost

The first term ||x[k-N] - x_bar||^2_{P_bar^{-1}} is the *arrival cost*. It
summarises all information from measurements before the estimation window. As
the window slides forward, the arrival cost is updated to preserve historical
information.

For linear systems with Gaussian noise, the exact arrival cost is the Kalman
filter covariance. For nonlinear systems, an EKF approximation is commonly
used -- ctrlpp auto-constructs a companion EKF for this purpose.

### Relationship to Kalman Filter

For linear systems without constraints, MHE with infinite window length
recovers the Kalman filter exactly. With finite window length, MHE
approximates the Kalman filter but can additionally enforce state constraints.

### QP/NLP Formulation

- **Linear MHE**: the dynamics are linear, yielding a quadratic program.
  ctrlpp uses OSQP as the default QP solver.
- **Nonlinear MHE**: nonlinear dynamics yield a nonlinear program. ctrlpp
  uses NLopt as the default NLP solver.

Both formulations support box constraints on states.

## References

- **Rao, C. V., Rawlings, J. B., and Mayne, D. Q.** "Constrained State
  Estimation for Nonlinear Discrete-Time Systems: Stability and Moving
  Horizon Approximations." *IEEE Trans. Automatic Control*, 48(2):246--258,
  2003. DOI: 10.1109/TAC.2002.808470.
  Establishes the theoretical foundation for MHE: stability, convergence, and
  the role of the arrival cost.

- **Diehl, M., Ferreau, H. J., and Haverbeke, N.** "Efficient Numerical
  Methods for Nonlinear MPC and Moving Horizon Estimation." In *Nonlinear
  Model Predictive Control*, LNCIS vol. 384, pp. 391--417. Springer, 2009.
  DOI: 10.1007/978-3-642-01094-1_32.
  Covers computational aspects of MHE including real-time iteration schemes.

## Related API Pages

- [mhe](../mpc/mhe.md) -- linear moving horizon estimation (OSQP)
- [nmhe](../mpc/nmhe.md) -- nonlinear moving horizon estimation (NLopt)
- [MPC Theory](mpc-theory.md) -- the control dual of MHE
