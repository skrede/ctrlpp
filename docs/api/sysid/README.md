# System Identification

Online and offline methods for identifying dynamic system models from
input-output data. Online methods (RLS, recursive ARX) update estimates
sample-by-sample; offline methods (batch ARX, N4SID) operate on collected
datasets.

## Types

### Online

- [rls](rls.md) -- Recursive least squares with bounded covariance
- [recursive_arx](recursive-arx.md) -- Recursive ARX with `to_state_space()` conversion

### Offline

- [batch_arx](batch-arx.md) -- Batch ARX identification via QR factorisation
- [n4sid](n4sid.md) -- Subspace identification via BDCSVD

### Utilities

- [fit_metrics](fit-metrics.md) -- NRMSE, VAF, and other goodness-of-fit metrics
- [sysid_result](sysid-result.md) -- Identification result container

## When to use

Pick **RLS** for scalar-output online parameter tracking with exponential
forgetting.

Pick **recursive ARX** when you need an online ARX model that can be converted
to state-space for direct use with Kalman filters or MPC.

Pick **batch ARX** for offline identification of ARX models from a collected
dataset -- uses QR factorisation for numerical stability.

Pick **N4SID** for offline subspace identification when you want a state-space
model directly without specifying model orders.

## Theory

For the mathematical background, see
[System Identification Theory](../reference/sysid-theory.md).
