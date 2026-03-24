# batch_arx

Batch ARX model identification via QR decomposition. Given input/output data sequences, estimates an ARX(NA, NB) model and returns the identified discrete state-space system in observer canonical form along with fit metrics (NRMSE and VAF).

## Header and Alias

| Form | Header |
|------|--------|
| `batch_arx<NA, NB>(Y, U)` | `#include <ctrlpp/sysid/batch_arx.h>` |
| (convenience) | `#include <ctrlpp/sysid.h>` |

```cpp
template <std::size_t NA, std::size_t NB, typename Derived1, typename Derived2>
arx_result<typename Derived1::Scalar, NA, 1, 1>
batch_arx(const Eigen::MatrixBase<Derived1>& Y,
          const Eigen::MatrixBase<Derived2>& U);
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

## Return Type

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct arx_result {
    discrete_state_space<Scalar, NX, NU, NY> system;
    fit_metrics<Scalar> metrics;
};
```

The identified system is in observer canonical form with `NX = NA` states. The `metrics` field contains NRMSE and VAF computed by simulating the identified model against the original data.

## Usage Example

```cpp
#include "ctrlpp/sysid/batch_arx.h"

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

    std::cout << "Identified system:\n"
              << "  A = " << result.system.A << "\n"
              << "  B = " << result.system.B << "\n"
              << "  NRMSE = " << result.metrics.nrmse << "\n"
              << "  VAF   = " << result.metrics.vaf << " %\n";
}
```

## See Also

- [recursive-arx](recursive-arx.md) -- online recursive variant
- [n4sid](n4sid.md) -- subspace identification
- [fit-metrics](fit-metrics.md) -- goodness-of-fit metrics
- [sysid-result](sysid-result.md) -- result container types
