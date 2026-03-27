# System Identification

System identification builds mathematical models of dynamic systems from
measured input-output data. The goal is to find a model that accurately
predicts the system's behaviour, suitable for control design (LQR, MPC) or
state estimation (Kalman filter, MHE)
[1, Ch. 1, pp. 1--18].

System identification is essential when first-principles models are
unavailable, too complex, or insufficiently accurate. It provides a
data-driven alternative that can capture the actual behaviour of physical
systems including unmodelled dynamics and nonlinearities.

## The Identification Problem

Given a sequence of input-output measurements $\{u(t), y(t)\}_{t=1}^{N}$,
the identification problem is to find a model $\hat{y}(t) = g(u(t), \theta)$
parameterised by $\theta$ that minimises the prediction error
[1, Ch. 7, pp. 197--234]:

$$
V(\theta) = \sum_{t=1}^{N} \lVert y(t) - \hat{y}(t \mid \theta) \rVert^2
$$

The model structure (how $\hat{y}$ depends on past data and $\theta$) and the
model order (how many parameters) must be chosen by the user.

## ARX Models

The AutoRegressive with eXogenous input (ARX) model relates the current output
to past outputs and inputs [1, Ch. 4, pp. 83--108]:

$$
y(t) = a_1 y(t-1) + \cdots + a_{n_a} y(t-n_a) + b_1 u(t-1) + \cdots + b_{n_b} u(t-n_b) + e(t)
$$

In compact regression form with parameter vector
$\theta = [a_1, \ldots, a_{n_a}, b_1, \ldots, b_{n_b}]^\top$ and regressor
$\varphi(t) = [y(t-1), \ldots, y(t-n_a), u(t-1), \ldots, u(t-n_b)]^\top$
[1, eq. (4.2), p. 86]:

$$
y(t) = \varphi(t)^\top \theta + e(t)
$$

ARX is linear-in-parameters: the unknown coefficients can be estimated by
ordinary least squares. The model orders $n_a$ (output) and $n_b$ (input)
are selected by the user, typically using cross-validation or information
criteria (AIC, BIC) [1, Ch. 16, pp. 487--514].

## Batch Least Squares

For offline identification with $N$ data samples, the ARX parameters are
estimated by minimising the sum of squared residuals
[1, Ch. 7, pp. 199--205]:

$$
\hat{\theta} = \arg\min_\theta \sum_{t=1}^{N} (y(t) - \varphi(t)^\top \theta)^2
$$

The solution is the well-known least squares formula:

$$
\hat{\theta} = (\Phi^\top \Phi)^{-1} \Phi^\top Y
$$

where $\Phi = [\varphi(1), \ldots, \varphi(N)]^\top$ is the regression
matrix and $Y = [y(1), \ldots, y(N)]^\top$. For numerical stability, QR
factorisation is preferred over the normal equations
[2, Sec. 6.3, pp. 203--208].

## Recursive Least Squares (RLS)

RLS updates parameter estimates online as new data arrives, using the
Sherman-Morrison rank-1 covariance update
[1, Ch. 11, pp. 361--386]:

$$
K_t = \frac{P_{t-1} \, \varphi_t}{\lambda + \varphi_t^\top P_{t-1} \, \varphi_t}
$$

$$
\hat{\theta}_t = \hat{\theta}_{t-1} + K_t \bigl(y_t - \varphi_t^\top \hat{\theta}_{t-1}\bigr)
$$

$$
P_t = \frac{1}{\lambda} \bigl(P_{t-1} - K_t \, \varphi_t^\top P_{t-1}\bigr)
$$

### Forgetting Factor

The forgetting factor $\lambda \in (0, 1]$ discounts old data, allowing the
estimator to track slowly varying parameters
[1, Sec. 11.3, pp. 371--376]:

- $\lambda = 1$: all data weighted equally (standard least squares)
- $\lambda < 1$: exponential discounting of old data
- Typical values: $\lambda = 0.95$--$0.99$

Smaller $\lambda$ increases sensitivity to recent data but reduces estimation
precision and can cause covariance blow-up. Bounding the covariance matrix $P$
(e.g., resetting when its trace exceeds a threshold) prevents numerical
instability [1, Sec. 11.3, p. 374].

## Subspace Identification (N4SID)

Subspace identification estimates a state-space model $(A, B, C, D)$ directly
from input-output data without requiring explicit model order specification.
The N4SID algorithm [3, pp. 75--93]:

### Block Hankel Matrices

1. Construct block Hankel matrices from the data. For output sequence
   $\{y_0, y_1, \ldots\}$ with block row dimension $i$:

$$
Y = \begin{bmatrix} y_0 & y_1 & \cdots & y_{j-1} \\ y_1 & y_2 & \cdots & y_j \\ \vdots & \vdots & \ddots & \vdots \\ y_{i-1} & y_i & \cdots & y_{i+j-2} \end{bmatrix}
$$

### Observability Matrix Estimation

2. Estimate the column space of the extended observability matrix
   $\mathcal{O}_i = [C^\top, (CA)^\top, \ldots, (CA^{i-1})^\top]^\top$ via
   singular value decomposition (SVD). The model order $n$ is determined by
   the gap in the singular value spectrum [3, Sec. 3, pp. 80--85].

### State-Space Extraction

3. Extract $(A, B, C, D)$ from the identified subspace using least-squares.
   The matrices are unique up to a similarity transformation
   [3, Sec. 4, pp. 85--90].

N4SID handles multi-input multi-output (MIMO) systems naturally and
automatically determines the model order from the data.

## ARX to State-Space Conversion

Identified ARX parameters can be converted to a discrete state-space
representation in observer canonical form
[1, Ch. 4, pp. 95--100]:

$$
A = \begin{bmatrix} -a_1 & 1 & 0 & \cdots & 0 \\ -a_2 & 0 & 1 & \cdots & 0 \\ \vdots & & & \ddots & \vdots \\ -a_{n_a-1} & 0 & 0 & \cdots & 1 \\ -a_{n_a} & 0 & 0 & \cdots & 0 \end{bmatrix}, \quad
B = \begin{bmatrix} b_1 \\ b_2 \\ \vdots \\ b_{n_a} \end{bmatrix}, \quad
C = \begin{bmatrix} 1 & 0 & \cdots & 0 \end{bmatrix}
$$

This state-space form can be directly used with Kalman filters, LQR, and MPC.

## Model Validation

Fit metrics quantify how well the identified model matches held-out data
[1, Ch. 16, pp. 492--498]:

- **NRMSE** (Normalised Root Mean Square Error):

$$
\text{NRMSE} = 1 - \frac{\lVert y - \hat{y} \rVert}{\lVert y - \bar{y} \rVert}
$$

  A value of 1.0 is a perfect fit; values above 0.8 are typically good.

- **VAF** (Variance Accounted For):

$$
\text{VAF} = 100 \left(1 - \frac{\text{var}(y - \hat{y})}{\text{var}(y)}\right) \%
$$

  Measures the percentage of output variance explained by the model.

Always validate on data not used for identification (cross-validation) to
detect overfitting.

## References

[1] L. Ljung, "System Identification: Theory for the User," 2nd ed., Prentice
Hall, 1999.

[2] A. K. Tangirala, "Principles of System Identification: Theory and
Practice," 2nd ed., CRC Press, 2018.

[3] P. Van Overschee and B. De Moor, "N4SID: Subspace Algorithms for the
Identification of Combined Deterministic-Stochastic Systems," Automatica,
vol. 30, no. 1, pp. 75--93, 1994.
