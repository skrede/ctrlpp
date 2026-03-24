# System Identification Theory

System identification builds mathematical models from measured input-output
data. The goal is to find a model that accurately predicts the system's
behaviour, suitable for control design (LQR, MPC) or state estimation
(Kalman filter).

## Key Concepts

### ARX Models

The AutoRegressive with eXogenous input (ARX) model relates current output to
past outputs and inputs:

```
y(t) = a_1 y(t-1) + ... + a_na y(t-na)
     + b_1 u(t-1) + ... + b_nb u(t-nb) + e(t)
```

ARX is a linear-in-parameters model: the unknown coefficients [a_1...a_na,
b_1...b_nb] can be estimated by ordinary least squares. The model orders
na and nb are chosen by the user.

Once identified, the ARX parameters are converted to a discrete state-space
model in observer canonical form for use with ctrlpp's controllers.

### Subspace Methods (N4SID)

Subspace identification estimates a state-space model directly from
input-output data without requiring model order specification. The N4SID
algorithm:

1. Constructs block Hankel matrices from the data
2. Estimates the column space of the observability matrix via SVD
3. Extracts (A, B, C, D) from the identified subspace

The model order is determined by the singular value gap. N4SID handles
multi-input multi-output (MIMO) systems naturally.

### Recursive Identification (RLS)

Recursive Least Squares (RLS) updates parameter estimates online as new data
arrives:

```
K[t]     = P[t-1] phi[t] / (lambda + phi[t]^T P[t-1] phi[t])
theta[t] = theta[t-1] + K[t] (y[t] - phi[t]^T theta[t-1])
P[t]     = (P[t-1] - K[t] phi[t]^T P[t-1]) / lambda
```

The forgetting factor lambda (0 < lambda <= 1) discounts old data, allowing
the estimator to track slowly varying parameters. ctrlpp uses the
Sherman-Morrison formula for efficient rank-1 covariance updates.

### Model Validation

Fit metrics quantify how well the identified model matches the data:

- **NRMSE** (Normalised Root Mean Square Error): 1 - ||y - y_hat|| / ||y - mean(y)||.
  A value of 1.0 is a perfect fit; values above 0.8 are typically good.
- **VAF** (Variance Accounted For): measures the percentage of output variance
  explained by the model. 100% is perfect.

Always validate on data not used for identification (cross-validation) to
check for overfitting.

## References

- **Ljung, L.** *System Identification: Theory for the User.* Prentice Hall,
  2nd ed., 1999. ISBN 978-0-13-656695-3.
  The standard reference for system identification. Covers ARX, ARMAX, output
  error, subspace methods, and recursive estimation.

- **Van Overschee, P. and De Moor, B.** "N4SID: Subspace Algorithms for the
  Identification of Combined Deterministic-Stochastic Systems." *Automatica*,
  30(1):75--93, 1994. DOI: 10.1016/0005-1098(94)90230-5.
  The original N4SID paper establishing subspace identification via oblique
  projections.

## Related API Pages

- [batch_arx](../sysid/batch-arx.md) -- offline batch ARX identification (QR)
- [n4sid](../sysid/n4sid.md) -- subspace identification (BDCSVD)
- [rls](../sysid/rls.md) -- recursive least squares
- [recursive_arx](../sysid/recursive-arx.md) -- recursive ARX identification
- [fit_metrics](../sysid/fit-metrics.md) -- NRMSE and VAF
- [Sysid Workflow Guide](../guides/sysid/workflow.md) -- identify, model,
  and control pipeline
