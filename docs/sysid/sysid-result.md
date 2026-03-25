# sysid_result

Result containers returned by the system identification algorithms. Each result type bundles the identified state-space model with algorithm-specific metadata and fit metrics.

## Header and Alias

| Form | Header |
|------|--------|
| `arx_result<Scalar, NX, NU, NY>` | `#include <ctrlpp/sysid/sysid_result.h>` |
| `n4sid_result<Scalar, NX, NU, NY>` | `#include <ctrlpp/sysid/sysid_result.h>` |
| (convenience) | `#include <ctrlpp/sysid.h>` |

## arx_result

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct arx_result {
    discrete_state_space<Scalar, NX, NU, NY> system;
    fit_metrics<Scalar> metrics;
};
```

Returned by `batch_arx`. Contains the identified system in observer canonical form and fit metrics computed by simulating the model against the training data.

| Field | Type | Description |
|-------|------|-------------|
| `system` | `discrete_state_space<Scalar, NX, NU, NY>` | Identified state-space model |
| `metrics` | `fit_metrics<Scalar>` | NRMSE and VAF |

## n4sid_result

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct n4sid_result {
    discrete_state_space<Scalar, NX, NU, NY> system;
    Eigen::VectorX<Scalar> singular_values;
    fit_metrics<Scalar> metrics;
    Scalar condition_number;
};
```

Returned by `n4sid`. Contains the identified system, the oblique projection singular values (for model order selection), fit metrics, and the condition number of the observability matrix truncation.

| Field | Type | Description |
|-------|------|-------------|
| `system` | `discrete_state_space<Scalar, NX, NU, NY>` | Identified state-space model |
| `singular_values` | `Eigen::VectorX<Scalar>` | SVD singular values from oblique projection |
| `metrics` | `fit_metrics<Scalar>` | NRMSE and VAF |
| `condition_number` | `Scalar` | Condition number of the truncated observability matrix |

## Usage Example

```cpp
#include "ctrlpp/sysid/batch_arx.h"
#include "ctrlpp/sysid/n4sid.h"

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
    std::cout << "ARX: NRMSE=" << arx.metrics.nrmse
              << "  VAF=" << arx.metrics.vaf << "%\n";

    // N4SID identification
    auto ss = ctrlpp::n4sid<1>(Y, U);
    std::cout << "N4SID: NRMSE=" << ss.metrics.nrmse
              << "  VAF=" << ss.metrics.vaf << "%"
              << "  cond=" << ss.condition_number << "\n";
}
```

## See Also

- [fit-metrics](fit-metrics.md) -- goodness-of-fit metric computation
- [batch-arx](batch-arx.md) -- batch ARX identification
- [n4sid](n4sid.md) -- subspace identification
- [recursive-arx](recursive-arx.md) -- recursive ARX identification
- [rls](rls.md) -- recursive least squares
- [model/state-space](../model/state-space.md) -- state-space representation
