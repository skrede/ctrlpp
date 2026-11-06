# System Identification

Online and offline methods for identifying linear dynamic system models from
input-output data. Online methods (RLS, recursive ARX) update estimates
sample-by-sample; offline methods (batch ARX, MOESP) operate on collected
datasets. The ARX and MOESP routines identify single-input single-output
systems at orders the caller chooses.

`#include <ctrlpp/sysid.h>` includes every public header of the module.

## Types

### Online

- [rls](rls.md)<br/> Recursive least squares with bounded covariance
- [recursive_arx](recursive-arx.md)<br/> Recursive ARX with `to_state_space()` conversion

### Offline

- [batch_arx](batch-arx.md)<br/> Batch ARX identification via column-pivoting QR factorization
- [moesp](moesp.md)<br/> PO-MOESP subspace identification at a caller-chosen order

### Utilities

- [fit_metrics](fit-metrics.md)<br/> NRMSE and VAF goodness-of-fit metrics
- [sysid_result](sysid-result.md)<br/> Identification result containers
- `sysid_error`, the record-validation error of `batch_arx` and `moesp`, is declared in `ctrlpp/sysid/sysid_types.h` and listed on [batch_arx](batch-arx.md#record-validation-sysid_error)

## When to use

Pick **RLS** for scalar-output online parameter tracking with exponential
forgetting.

Pick **recursive ARX** when you need an online ARX model that can be converted
to state-space for direct use with Kalman filters or MPC.

Pick **batch ARX** for offline identification of ARX models from a collected
dataset &mdash; uses QR factorization for numerical stability.

Pick **MOESP** for offline subspace identification when you want a state-space
model directly. You choose the order `NX`; `moesp_singular_values` returns the
singular-value spectrum to inform that choice.

## Theory

For the mathematical background, see
[System Identification Theory](../../background/sysid.md).
