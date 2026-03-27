# n4sid

Subspace system identification using the N4SID algorithm. Identifies a discrete-time linear state-space model of a given order from input/output data via oblique projection and singular value decomposition. Returns the identified system in state-space form with fit metrics and condition number.

## Header and Alias

| Form | Header |
|------|--------|
| `n4sid<NX>(Y, U, block_rows)` | `#include <ctrlpp/sysid/n4sid.h>` |
| (convenience) | `#include <ctrlpp/sysid.h>` |

```cpp
template <std::size_t NX, typename Derived1, typename Derived2>
n4sid_result<typename Derived1::Scalar, NX, 1, 1>
n4sid(const Eigen::MatrixBase<Derived1>& Y,
      const Eigen::MatrixBase<Derived2>& U,
      std::size_t block_rows = 0);
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `NX` | `>= 1` | Target state-space order |

## Function Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `Y` | `Eigen::MatrixBase` | -- | Output data as a 1-by-N row matrix |
| `U` | `Eigen::MatrixBase` | -- | Input data as a 1-by-N row matrix |
| `block_rows` | `std::size_t` | `0` | Block Hankel matrix row count. When 0, defaults to `min(N/4, 30)`. |

## Return Type

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct n4sid_result {
    discrete_state_space<Scalar, NX, NU, NY> system;
    Eigen::VectorX<Scalar> singular_values;
    fit_metrics<Scalar> metrics;
    Scalar condition_number;
};
```

The `singular_values` field from the oblique projection SVD can be used to determine the appropriate model order (look for a gap in the singular value spectrum).

## Singular Value Helper

```cpp
template <typename Derived1, typename Derived2>
Eigen::VectorX<typename Derived1::Scalar>
n4sid_singular_values(const Eigen::MatrixBase<Derived1>& Y,
                      const Eigen::MatrixBase<Derived2>& U,
                      std::size_t block_rows = 0);
```

Returns just the singular values without performing full identification. Useful for model order selection.

## Usage Example

```cpp
#include <ctrlpp/sysid/n4sid.h>
#include <ctrlpp/model/analysis.h>

#include <Eigen/Dense>

#include <iostream>
#include <random>

int main()
{
    // Generate data from a 2nd-order system
    constexpr int N = 1000;
    Eigen::RowVectorXd Y(N);
    Eigen::RowVectorXd U(N);

    Eigen::Matrix2d A;
    A << 0.9, 0.1, -0.2, 0.8;
    Eigen::Vector2d B(0.5, 0.3);
    Eigen::RowVector2d C(1.0, 0.0);

    std::mt19937 rng(42);
    std::normal_distribution<double> input_dist(0.0, 1.0);
    std::normal_distribution<double> noise(0.0, 0.01);

    Eigen::Vector2d x = Eigen::Vector2d::Zero();

    for(int t = 0; t < N; ++t)
    {
        U(t) = input_dist(rng);
        Y(t) = (C * x)(0) + noise(rng);
        Eigen::Matrix<double, 1, 1> u_vec;
        u_vec << U(t);
        x = A * x + B * U(t);
    }

    // Check singular values for order selection
    auto sv = ctrlpp::n4sid_singular_values(Y, U);
    std::cout << "Singular values: " << sv.transpose() << "\n\n";

    // Identify 2nd-order model
    auto result = ctrlpp::n4sid<2>(Y, U);

    std::cout << "Identified system:\n"
              << "  A =\n" << result.system.A << "\n"
              << "  B =\n" << result.system.B << "\n"
              << "  C = " << result.system.C << "\n"
              << "  NRMSE = " << result.metrics.nrmse << "\n"
              << "  VAF   = " << result.metrics.vaf << " %\n"
              << "  Cond  = " << result.condition_number << "\n";
}
```

## See Also

- [batch-arx](batch-arx.md) -- ARX identification
- [model/state-space](../model/state-space.md) -- state-space representation
- [sysid-result](sysid-result.md) -- result container types
- [fit-metrics](fit-metrics.md) -- NRMSE and VAF metrics
- [guides/sysid/workflow](../../guides/sysid/workflow.md) -- system identification workflow guide
- [background/sysid](../../background/sysid.md) -- sysid theory and background
