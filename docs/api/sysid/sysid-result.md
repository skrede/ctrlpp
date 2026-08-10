# sysid_result

Result containers returned by the system identification algorithms. Each result type bundles the identified state-space model with algorithm-specific metadata and fit metrics.

## Header and Alias

| Form | Header |
|------|--------|
| `arx_result<Scalar, NX, NU, NY>` | `#include <ctrlpp/sysid/sysid_result.h>` |
| `arx_diagnostics<Scalar>` | `#include <ctrlpp/sysid/sysid_result.h>` |
| `moesp_result<Scalar, NX, NU, NY>` | `#include <ctrlpp/sysid/sysid_result.h>` |
| (convenience) | `#include <ctrlpp/sysid.h>` |

## arx_result

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct arx_result {
    discrete_state_space<Scalar, NX, NU, NY> system;
    fit_metrics<Scalar> metrics;
    arx_diagnostics<Scalar> diagnostics;
};
```

Carried on the value branch of the `ctrlpp::expected<arx_result<...>, sysid_error>` returned by `batch_arx`. Contains the identified system in observer canonical form, fit metrics computed by simulating the model against the training data, and diagnostics from the least-squares stage that produced it. See [batch-arx](batch-arx.md) for the rejection list on the error branch.

| Field | Type | Description |
|-------|------|-------------|
| `system` | `discrete_state_space<Scalar, NX, NU, NY>` | Identified state-space model |
| `metrics` | `fit_metrics<Scalar>` | NRMSE and VAF |
| `diagnostics` | `arx_diagnostics<Scalar>` | Rank and residual of the least-squares fit |

## arx_diagnostics

```cpp
template <typename Scalar>
struct arx_diagnostics {
    std::size_t effective_samples;
    std::size_t parameter_count;
    std::size_t numerical_rank;
    Scalar residual_norm;
};
```

What the least-squares stage resolved, as distinct from how well the identified model reproduces the record. `metrics` answers "does this model predict the data"; `diagnostics` answers "did the data determine this model".

| Field | Type | Description |
|-------|------|-------------|
| `effective_samples` | `std::size_t` | Regressor rows the fit was formed from: the record length less `max(NA, NB)` |
| `parameter_count` | `std::size_t` | Regressor columns, `NA + NB` |
| `numerical_rank` | `std::size_t` | Rank of the regressor matrix, counted by the column-pivoting QR against the active threshold |
| `residual_norm` | `Scalar` | Euclidean norm of the least-squares residual over those rows |

A `numerical_rank` below `parameter_count` means the record did not determine every coefficient direction. Such a fit is returned rather than refused, so this field is what makes the deficiency visible; see [batch-arx](batch-arx.md) for the rank threshold that decides it.

`effective_samples` is smaller than the scored sample count. The goodness-of-fit metrics simulate the identified model over the whole record from a zero initial state, so they score the `max(NA, NB)` startup samples that the regressor skips.

`residual_norm` is the residual of the **fit**, over the same rows the rank was counted from, not the residual of the simulation the metrics score. On a record generated exactly from a model of the fitted order it sits at the arithmetic's rounding level; on a disturbed record it cannot exceed the norm of the disturbance, because the fit returns a minimizer.

## moesp_result

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct moesp_result {
    discrete_state_space<Scalar, NX, NU, NY> system;
    Eigen::VectorX<Scalar> singular_values;
    fit_metrics<Scalar> metrics;
    Scalar condition_number;
};
```

Returned by `moesp`. Contains the identified system, the oblique projection singular values (for model order selection), fit metrics, and the condition number of the observability matrix truncation.

| Field | Type | Description |
|-------|------|-------------|
| `system` | `discrete_state_space<Scalar, NX, NU, NY>` | Identified state-space model |
| `singular_values` | `Eigen::VectorX<Scalar>` | SVD singular values from oblique projection |
| `metrics` | `fit_metrics<Scalar>` | NRMSE and VAF |
| `condition_number` | `Scalar` | Condition number of the truncated observability matrix |

## Usage Example

```cpp
#include <ctrlpp/sysid/batch_arx.h>
#include <ctrlpp/sysid/moesp.h>

#include <Eigen/Dense>

#include <iostream>
#include <random>

int main()
{
    constexpr int N = 500;
    Eigen::RowVectorXd Y(N);
    Eigen::RowVectorXd U(N);

    std::mt19937 rng(42);
    std::normal_distribution<double> input_dist(0.0, 1.0);

    Y(0) = 0.0;
    U(0) = input_dist(rng);
    for(int t = 1; t < N; ++t)
    {
        U(t) = input_dist(rng);
        Y(t) = 0.8 * Y(t - 1) + 0.4 * U(t - 1);
    }

    // ARX identification
    auto arx = ctrlpp::batch_arx<1, 1>(Y, U);
    if(!arx)
    {
        std::cerr << "ARX identification rejected the record\n";
        return 1;
    }
    std::cout << "ARX: NRMSE=" << arx->metrics.nrmse
              << "  VAF=" << arx->metrics.vaf << "%\n";
    if(arx->diagnostics.numerical_rank < arx->diagnostics.parameter_count)
        std::cout << "  warning: the record determined only "
                  << arx->diagnostics.numerical_rank << " of "
                  << arx->diagnostics.parameter_count << " coefficients\n";

    // MOESP identification
    auto ss = ctrlpp::moesp<1>(Y, U);
    std::cout << "MOESP: NRMSE=" << ss.metrics.nrmse
              << "  VAF=" << ss.metrics.vaf << "%"
              << "  cond=" << ss.condition_number << "\n";
}
```

## See Also

- [fit-metrics](fit-metrics.md)<br/> goodness-of-fit metric computation
- [batch-arx](batch-arx.md)<br/> batch ARX identification
- [moesp](moesp.md)<br/> subspace identification
- [recursive-arx](recursive-arx.md)<br/> recursive ARX identification
- [rls](rls.md)<br/> recursive least squares
- [model/state-space](../model/state-space.md)<br/> state-space representation
