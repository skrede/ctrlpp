# moesp

Subspace system identification using the PO-MOESP algorithm. Identifies a discrete-time linear state-space model of a given order from input/output data by projecting the future outputs onto the past inputs and outputs (the instrumental variables) and the orthogonal complement of the future inputs, then recovering the extended observability matrix via a singular value decomposition. Returns the identified system in state-space form with fit metrics and condition number.

## Header and Alias

| Form | Header |
|------|--------|
| `moesp<NX>(Y, U, block_rows)` | `#include <ctrlpp/sysid/moesp.h>` |
| (convenience) | `#include <ctrlpp/sysid.h>` |

```cpp
template <std::size_t NX, typename Derived1, typename Derived2>
moesp_result<typename Derived1::Scalar, NX, 1, 1>
moesp(const Eigen::MatrixBase<Derived1>& Y,
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
| `Y` | `Eigen::MatrixBase` | &mdash; | Output data as a 1-by-N row matrix |
| `U` | `Eigen::MatrixBase` | &mdash; | Input data as a 1-by-N row matrix |
| `block_rows` | `std::size_t` | `0` | Block Hankel matrix row count. When 0, defaults to `min(N/4, 30)`. |

## Return Type

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct moesp_result {
    discrete_state_space<Scalar, NX, NU, NY> system;
    Eigen::VectorX<Scalar> singular_values;
    fit_metrics<Scalar> metrics;
    Scalar condition_number;
};
```

The `singular_values` field from the oblique projection SVD can be used to determine the appropriate model order (look for a gap in the singular value spectrum). When the input data is rank-deficient or produces non-finite system matrices, `condition_number` is set to infinity and the system matrices will be default-initialized (zero). Check `std::isinf(result.condition_number)` to detect identification failure.

## Singular Value Helper

```cpp
template <typename Derived1, typename Derived2>
Eigen::VectorX<typename Derived1::Scalar>
moesp_singular_values(const Eigen::MatrixBase<Derived1>& Y,
                      const Eigen::MatrixBase<Derived2>& U,
                      std::size_t block_rows = 0);
```

Returns just the singular values without performing full identification. Useful for model order selection.

## Usage Example

```cpp
#include <ctrlpp/sysid/moesp.h>
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
    auto sv = ctrlpp::moesp_singular_values(Y, U);
    std::cout << "Singular values: " << sv.transpose() << "\n\n";

    // Identify 2nd-order model
    auto result = ctrlpp::moesp<2>(Y, U);

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

- [batch-arx](batch-arx.md)<br/> ARX identification
- [model/state-space](../model/state-space.md)<br/> state-space representation
- [sysid-result](sysid-result.md)<br/> result container types
- [fit-metrics](fit-metrics.md)<br/> NRMSE and VAF metrics
- [guides/sysid/workflow](../../guides/sysid/workflow.md)<br/> system identification workflow guide
- [background/sysid](../../background/sysid.md)<br/> sysid theory and background
