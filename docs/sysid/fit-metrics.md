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
fit_metrics<typename DerivedA::Scalar>
compute_fit_metrics(const Eigen::MatrixBase<DerivedA>& y_actual,
                    const Eigen::MatrixBase<DerivedB>& y_predicted);
```

Computes both metrics from two Eigen column vectors of the same length.

**NRMSE:** `||y_actual - y_predicted|| / ||y_actual - mean(y_actual)||`. A value of 0 indicates a perfect fit; values above 1 indicate the model is worse than predicting the mean.

**VAF:** `(1 - var(error) / var(y_actual)) * 100`. A value of 100% indicates the model explains all variance; values near 0% indicate no explanatory power.

Edge cases: constant signals (zero variance) return NRMSE = 0 and VAF = 100% when the prediction is also perfect, or infinity/-infinity otherwise.

## Usage Example

```cpp
#include "ctrlpp/sysid/fit_metrics.h"

#include <Eigen/Dense>

#include <iostream>

int main()
{
    Eigen::VectorXd y_actual(5);
    y_actual << 1.0, 2.0, 3.0, 4.0, 5.0;

    Eigen::VectorXd y_predicted(5);
    y_predicted << 1.1, 1.9, 3.2, 3.8, 5.1;

    auto metrics = ctrlpp::compute_fit_metrics(y_actual, y_predicted);

    std::cout << "NRMSE = " << metrics.nrmse << "\n"
              << "VAF   = " << metrics.vaf << " %\n";
}
```

## See Also

- [sysid-result](sysid-result.md) -- result containers that include fit_metrics
- [batch-arx](batch-arx.md) -- batch ARX identification
- [n4sid](n4sid.md) -- subspace identification
- [rls](rls.md) -- recursive least squares
- [reference/sysid-theory](../reference/sysid-theory.md) -- sysid theory and background
