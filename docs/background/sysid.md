# System Identification Theory

System identification builds mathematical models from measured input-output
data. The goal is to find a model that accurately predicts the system's
behaviour, suitable for control design (LQR, MPC) or state estimation
(Kalman filter).

## ARX Models

The AutoRegressive with eXogenous input (ARX) model relates current output to
past outputs and inputs:

$$
y(t) = a_1 y(t-1) + \cdots + a_{n_a} y(t - n_a) + b_1 u(t-1) + \cdots + b_{n_b} u(t - n_b) + e(t)
$$

In compact regression form with parameter vector
$\theta = [a_1, \ldots, a_{n_a}, b_1, \ldots, b_{n_b}]^\top$ and regressor
$\varphi(t) = [y(t-1), \ldots, y(t-n_a), u(t-1), \ldots, u(t-n_b)]^\top$:

$$
y(t) = \varphi(t)^\top \theta + e(t)
$$

ARX is linear-in-parameters: the unknown coefficients can be estimated by
ordinary least squares. The model orders $n_a$ and $n_b$ are chosen by the
user.

Once identified, the ARX parameters are converted to a discrete state-space
model in observer canonical form for use with ctrlpp's controllers and
estimators.

## Subspace Methods (N4SID)

Subspace identification estimates a state-space model $(A, B, C, D)$ directly
from input-output data without requiring model order specification. The N4SID
algorithm:

1. Constructs block Hankel matrices from the data. For output sequence
   $\{y_0, y_1, \ldots\}$ with block row dimension $i$:

$$
Y = \begin{bmatrix} y_0 & y_1 & \cdots & y_{j-1} \\ y_1 & y_2 & \cdots & y_j \\ \vdots & \vdots & \ddots & \vdots \\ y_{i-1} & y_i & \cdots & y_{i+j-2} \end{bmatrix}
$$

2. Estimates the column space of the extended observability matrix
   $\mathcal{O}_i = [C^\top, (CA)^\top, \ldots, (CA^{i-1})^\top]^\top$ via
   singular value decomposition (SVD). The model order $n$ is determined by
   the singular value gap.

3. Extracts $(A, B, C, D)$ from the identified subspace using least-squares.

N4SID handles multi-input multi-output (MIMO) systems naturally. ctrlpp uses
the BDCSVD algorithm for efficient computation.

## Recursive Least Squares (RLS)

RLS updates parameter estimates online as new data arrives, using the
Sherman-Morrison rank-1 covariance update:

$$
K_t = \frac{P_{t-1} \, \varphi_t}{\lambda + \varphi_t^\top P_{t-1} \, \varphi_t}
$$

$$
\hat{\theta}_t = \hat{\theta}_{t-1} + K_t \bigl(y_t - \varphi_t^\top \hat{\theta}_{t-1}\bigr)
$$

$$
P_t = \frac{1}{\lambda} \bigl(P_{t-1} - K_t \, \varphi_t^\top P_{t-1}\bigr)
$$

The forgetting factor $\lambda \in (0, 1]$ discounts old data, allowing the
estimator to track slowly varying parameters. When $\lambda = 1$, all data
is weighted equally (standard least squares). Smaller $\lambda$ increases
sensitivity to recent data but reduces estimation precision.

ctrlpp bounds the covariance matrix $P$ to prevent numerical blow-up when
$\lambda < 1$.

## Model Validation

Fit metrics quantify how well the identified model matches the data:

- **NRMSE** (Normalised Root Mean Square Error):

$$
\text{NRMSE} = 1 - \frac{\lVert y - \hat{y} \rVert}{\lVert y - \bar{y} \rVert}
$$

  A value of 1.0 is a perfect fit; values above 0.8 are typically good.

- **VAF** (Variance Accounted For):

$$
\text{VAF} = 100 \left(1 - \frac{\text{var}(y - \hat{y})}{\text{var}(y)}\right) \%
$$

  Measures the percentage of output variance explained by the model. 100% is
  perfect.

Always validate on data not used for identification (cross-validation) to
check for overfitting.

## References

- Ljung, L. (1999). *System Identification: Theory for the User.* Prentice
  Hall, 2nd ed. [`ljung1999`]
  The standard reference for system identification. Covers ARX, ARMAX,
  output error, subspace methods, and recursive estimation.

- Van Overschee, P. and De Moor, B. (1994). "N4SID: Subspace Algorithms for
  the Identification of Combined Deterministic-Stochastic Systems."
  *Automatica*, 30(1):75--93. [`vanoverschee1994`]
  The original N4SID paper establishing subspace identification via oblique
  projections.

## Related API Pages

- [batch_arx](../api/sysid/batch-arx.md) -- offline batch ARX identification (QR)
- [n4sid](../api/sysid/n4sid.md) -- subspace identification (BDCSVD)
- [rls](../api/sysid/rls.md) -- recursive least squares
- [recursive_arx](../api/sysid/recursive-arx.md) -- recursive ARX identification
- [fit_metrics](../api/sysid/fit-metrics.md) -- NRMSE and VAF
- [Sysid Workflow Guide](../guides/sysid/workflow.md) -- identify, model,
  and control pipeline
