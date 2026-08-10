# batch_arx

Batch ARX model identification via QR decomposition. Given input/output data sequences, estimates an ARX(NA, NB) model and returns the identified discrete state-space system in observer canonical form along with fit metrics (NRMSE and VAF) and diagnostics from the least-squares stage.

The routine validates its data records and returns `ctrlpp::expected<arx_result<...>, sysid_error>`: the identified model on success, or the specific `sysid_error` enumerator describing the rejected record. See [Record Validation: sysid_error](#record-validation-sysid_error).

## Header and Alias

| Form | Header |
|------|--------|
| `batch_arx<NA, NB>(Y, U)` | `#include <ctrlpp/sysid/batch_arx.h>` |
| (convenience) | `#include <ctrlpp/sysid.h>` |

```cpp
template <std::size_t NA, std::size_t NB, typename Derived1, typename Derived2>
auto
batch_arx(const Eigen::MatrixBase<Derived1>& Y,
          const Eigen::MatrixBase<Derived2>& U,
          std::optional<typename Derived1::Scalar> rank_threshold = std::nullopt)
    -> ctrlpp::expected<arx_result<typename Derived1::Scalar, std::max(NA, NB), 1, 1>,
                        sysid_error>;
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `NA` | `>= 1` | Number of auto-regressive (output) terms |
| `NB` | `>= 1` | Number of exogenous (input) terms |

## Function Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `Y` | `Eigen::MatrixBase` | Output data as a 1-by-N row matrix |
| `U` | `Eigen::MatrixBase` | Input data as a 1-by-N row matrix |
| `rank_threshold` | `std::optional<Scalar>` | Relative cutoff for the rank count. Omit for Eigen's own default |

## Return Type

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct arx_result {
    discrete_state_space<Scalar, NX, NU, NY> system;
    fit_metrics<Scalar> metrics;
    arx_diagnostics<Scalar> diagnostics;
};
```

The `arx_result` above is carried on the value branch of `ctrlpp::expected<arx_result<...>, sysid_error>`; the error branch carries a `sysid_error`.

The identified system is in observer canonical form with `NX = max(NA, NB)` states. The realization dimension is the larger of the auto-regressive and exogenous orders, so every b-coefficient is represented even when `NB > NA` (Ljung 1999, Ch. 4). The `metrics` field contains NRMSE and VAF computed by simulating the identified model against the original data, and the `diagnostics` field reports what the least-squares stage resolved; see [sysid-result](sysid-result.md).

## What the Fit Metrics Score

The metrics simulate the identified model **from a zero initial state** and score every sample of the record. The first `max(NA, NB)` predictions are therefore produced from a state the identification never estimated, and they are included in NRMSE and VAF rather than skipped. The regressor does skip those samples, so the row count the fit was actually formed from is smaller than the scored sample count; it is reported as `diagnostics.effective_samples`.

## Rank Deficiency and the Rank Threshold

A record that passes all four validity checks can still leave coefficient directions undetermined. Such a fit is **returned with the deficiency reported**, not refused: `diagnostics.numerical_rank` below `diagnostics.parameter_count` is the signal, and no `sysid_error` enumerator corresponds to it. The returned coefficients are finite in that case: the solve runs through the resolved pivots only and zeroes the rest.

`rank_threshold` is a **reporting policy and nothing else**. It is forwarded to the column-pivoting QR's rank test, which counts a pivot as resolved when it is strictly greater than this multiple of the largest pivot, so the value is relative and dimensionless. How many pivots the solve runs through is fixed during the factorization from machine epsilon alone, so no value of `rank_threshold` changes the identified system, the fit metrics or the residual norm. Passing `1.0` reports a rank of zero for a perfectly well-conditioned record and still returns that record's fit unchanged.

Omitting the argument reproduces Eigen's own default exactly: machine epsilon times the runtime diagonal size, which Eigen attributes to Higham's LDLT formula (`Eigen/src/QR/ColPivHouseholderQR.h`). That default counts the diagonal size only, so it does not grow with the record length, while the rounding a Householder QR leaves on a dependent pivot does. On a tall regressor an exactly collinear column can therefore still be counted as resolved. Where that matters, pass the backward-stable tolerance for a QR instead, the row count times machine epsilon:

```cpp
// The regressor row count: the record length less max(NA, NB), here 2.
auto const rows = static_cast<double>(Y.cols() - 2);
auto const result = ctrlpp::batch_arx<2, 2>(Y, U, rows * std::numeric_limits<double>::epsilon());
```

## Record Validation: sysid_error

Header: `#include <ctrlpp/sysid/sysid_types.h>` (pulled in by `batch_arx.h`)

```cpp
enum class sysid_error {
    record_length_mismatch,
    record_not_single_row,
    too_few_samples,
    non_finite_sample,
};
```

The checks run in the order listed, all of them before any index arithmetic, so a record with several defects reports the first matching enumerator:

| Enumerator | Rejected record |
|------------|-----------------|
| `record_length_mismatch` | `Y.cols() != U.cols()`. Sample `k` of one record is paired with sample `k` of the other, so unequal lengths have no consistent pairing |
| `record_not_single_row` | `Y.rows() != 1` or `U.rows() != 1`. The routine is single-input single-output and reads only row zero, so a multi-row record would be silently identified from a fraction of its data |
| `too_few_samples` | `Y.cols() <= max(NA, NB)`. One regressor row is formed per sample beyond `max(NA, NB)`, so a strictly greater sample count is exactly the condition for a nonempty regressor matrix. At the order the matrix is empty; below it the row count would be the result of an unsigned subtraction that wraps |
| `non_finite_sample` | A sample in either record is NaN or infinite. Samples enter the regressor and the least-squares solve directly, which propagates the value into every identified coefficient |

Every bound is an exact precondition of the routine's index arithmetic or data layout, not a tolerance. None of them is a quality judgement about the fit: `batch_arx` does **not** require the regressor row count to reach the parameter count `NA + NB`, so a record that passes all four checks can still yield an under-determined or rank-deficient fit. Judge that from the returned `arx_diagnostics`.

## Usage Example

```cpp
#include <ctrlpp/sysid/batch_arx.h>

#include <Eigen/Dense>

#include <iostream>
#include <random>

int main()
{
    // Generate data from a known system: y(t) = 0.7*y(t-1) + 0.3*u(t-1)
    constexpr int N = 500;
    Eigen::RowVectorXd Y(N);
    Eigen::RowVectorXd U(N);

    std::mt19937 rng(42);
    std::normal_distribution<double> input_dist(0.0, 1.0);
    std::normal_distribution<double> noise(0.0, 0.01);

    Y(0) = 0.0;
    U(0) = input_dist(rng);

    for(int t = 1; t < N; ++t)
    {
        U(t) = input_dist(rng);
        Y(t) = 0.7 * Y(t - 1) + 0.3 * U(t - 1) + noise(rng);
    }

    // Identify ARX(1, 1) model
    auto result = ctrlpp::batch_arx<1, 1>(Y, U);
    if(!result)
    {
        std::cerr << "identification rejected the record: "
                  << static_cast<int>(result.error()) << "\n";
        return 1;
    }

    if(result->diagnostics.numerical_rank < result->diagnostics.parameter_count)
    {
        std::cerr << "the record determined only "
                  << result->diagnostics.numerical_rank << " of "
                  << result->diagnostics.parameter_count << " coefficients\n";
        return 1;
    }

    std::cout << "Identified system:\n"
              << "  A = " << result->system.A << "\n"
              << "  B = " << result->system.B << "\n"
              << "  NRMSE = " << result->metrics.nrmse << "\n"
              << "  VAF   = " << result->metrics.vaf << " %\n"
              << "  fit residual = " << result->diagnostics.residual_norm
              << " over " << result->diagnostics.effective_samples << " rows\n";
}
```

## See Also

- [recursive-arx](recursive-arx.md)<br/> online recursive variant
- [moesp](moesp.md)<br/> subspace identification
- [fit-metrics](fit-metrics.md)<br/> goodness-of-fit metrics
- [sysid-result](sysid-result.md)<br/> result container types
- [guides/sysid/workflow](../../guides/sysid/workflow.md)<br/> system identification workflow guide
- [background/sysid](../../background/sysid.md)<br/> sysid theory and background
