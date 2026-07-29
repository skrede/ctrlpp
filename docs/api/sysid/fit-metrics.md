# fit_metrics

Goodness-of-fit metrics for evaluating identified system models. Computes NRMSE (Normalized Root Mean Square Error) and VAF (Variance Accounted For) from actual and predicted output sequences.

## Header and Alias

| Form | Header |
|------|--------|
| `fit_metrics<Scalar>` | `#include <ctrlpp/sysid/fit_metrics.h>` |
| (convenience) | `#include <ctrlpp/sysid.h>` |

## fit_metrics struct

```cpp
template <typename Scalar>
struct fit_metrics {
    Scalar nrmse;  // Normalized root mean square error (0 = perfect fit)
    Scalar vaf;    // Variance accounted for, in percent (100 = perfect fit)
};
```

## compute_fit_metrics

```cpp
template <typename DerivedA, typename DerivedB>
ctrlpp::expected<fit_metrics<typename DerivedA::Scalar>, fit_metrics_error>
compute_fit_metrics(const Eigen::MatrixBase<DerivedA>& y_actual,
                    const Eigen::MatrixBase<DerivedB>& y_predicted);
```

Computes both metrics from two Eigen column vectors of the same length.
Empty, mismatched, non-column, and non-finite records are rejected before any
reduction is evaluated.

**NRMSE:** `||y_actual - y_predicted|| / ||y_actual - mean(y_actual)||`. A value of 0 indicates a perfect fit; values above 1 indicate the model is worse than predicting the mean.

**VAF:** `(1 - var(error) / var(y_actual)) * 100`. A value of 100% indicates the model explains all variance; values near 0% indicate no explanatory power.

Edge cases: constant signals (zero variance) return NRMSE = 0 and VAF = 100% when the prediction is also perfect, or infinity/-infinity otherwise. Whether a quantity counts as vanishing is decided against a resolution floor **derived from the record's own magnitude** -- the counted rounding of the mean and the centring, `sqrt(n) * (n + 1)` operations, times the machine epsilon, times the largest sample -- and never against an absolute constant. An absolute one is wrong in both directions: at a record scale of 1e-17 it calls a genuinely varying record constant and reports a perfect fit for a predictor that explains nothing, and at a large scale it calls a constant record varying and then divides one rounding-level quantity by another.

## Usage Example

```cpp
#include <ctrlpp/sysid/fit_metrics.h>

#include <Eigen/Dense>

#include <iostream>

int main()
{
    Eigen::VectorXd y_actual(5);
    y_actual << 1.0, 2.0, 3.0, 4.0, 5.0;

    Eigen::VectorXd y_predicted(5);
    y_predicted << 1.1, 1.9, 3.2, 3.8, 5.1;

    auto metrics = ctrlpp::compute_fit_metrics(y_actual, y_predicted);
    if (!metrics) {
        return 1;
    }

    std::cout << "NRMSE = " << metrics->nrmse << "\n"
              << "VAF   = " << metrics->vaf << " %\n";
}
```

## See Also

- [sysid-result](sysid-result.md)<br/> result containers that include fit_metrics
- [batch-arx](batch-arx.md)<br/> batch ARX identification
- [moesp](moesp.md)<br/> subspace identification
- [rls](rls.md)<br/> recursive least squares
- [background/sysid](../../background/sysid.md)<br/> sysid theory and background
