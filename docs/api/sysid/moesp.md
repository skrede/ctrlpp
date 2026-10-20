# moesp

Subspace system identification using the PO-MOESP algorithm. Identifies a discrete-time linear state-space model of a given order from input/output data by projecting the future outputs onto the past inputs and outputs (the instrumental variables) and the orthogonal complement of the future inputs, then recovering the extended observability matrix via a singular value decomposition. `A` and `C` follow from the observability matrix; `B`, `D` and the initial state are then estimated together by one least-squares fit, so the record need not start at rest. Returns a `ctrlpp::expected` holding the identified system in state-space form with its estimated initial state, fit metrics and condition number, or a `sysid_error` naming why the record was refused.

## Header and Alias

| Form | Header |
|------|--------|
| `moesp<NX>(Y, U, block_rows)` | `#include <ctrlpp/sysid/moesp.h>` |
| (convenience) | `#include <ctrlpp/sysid.h>` |

```cpp
template <std::size_t NX, typename Derived1, typename Derived2>
auto moesp(const Eigen::MatrixBase<Derived1>& Y,
           const Eigen::MatrixBase<Derived2>& U,
           std::size_t block_rows = 0)
    -> ctrlpp::expected<moesp_result<typename Derived1::Scalar, NX, 1, 1>, sysid_error>;
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
    Vector<Scalar, NX> initial_state;
};
```

The `singular_values` field from the oblique projection SVD can be used to determine the appropriate model order (look for a gap in the singular value spectrum). `condition_number` is `sv(0) / sv(NX - 1)`, the ratio of the largest singular value to the last one the requested order consumes. It is a diagnostic of how well the order-NX subspace is separated, and on a returned model it is always finite.

`initial_state` is the state at the first sample, expressed in the identified state basis (the realization is determined only up to a similarity transform, so compare `C * initial_state` with the record's free response rather than with a state vector from another basis). It is estimated jointly with `B` and `D`: row `t` of the regression holds the responses `C X_t` of the `B` basis, the input `u_t`, and the free response `C A^t`, so a record that starts away from rest does not bias `B` and `D`. The fit `metrics` score the identified model simulated from `initial_state`, not from zero.

## Refusals

With `N` the sample count and `i` the block-row count actually used (`block_rows`, or `min(N/4, 30)` when it is zero), the record is checked in this order before anything is allocated. A record that violates two conditions reports the earlier one.

| Condition | Error | Why |
|-----------|-------|-----|
| `Y` and `U` hold different sample counts | `sysid_error::record_length_mismatch` | Sample k of one is paired with sample k of the other. |
| either record has a row count other than one | `sysid_error::record_not_single_row` | The routine is single-input single-output and reads only row zero. |
| `(i - 1) * ny < NX` | `sysid_error::block_rows_below_order` | The shift relation that recovers `A` has `(i - 1) * ny` rows and `NX` unknowns per column; below `NX` rows its minimum-norm solution is finite and wrong. The default block-row count is checked too, so a short record can report this without an explicit `block_rows`. |
| `N - 2i + 1 < i + NX` | `sysid_error::too_few_samples` | The block Hankel matrices cannot supply the `NX` projection columns the order needs. At `N - 2i + 1 <= 0` they do not exist at all. |
| a sample is NaN or infinite | `sysid_error::non_finite_sample` | Samples enter the factorization directly. |

After the projection, the identification is refused with `sysid_error::model_not_identified` when the condition number is not finite (the `NX`-th singular value is zero, as for an all-zero record, or so small against the leading one that the ratio overflows) or when an extracted system matrix or the estimated initial state is non-finite. No partially identified or zeroed model is returned.

## Singular Value Helper

```cpp
template <typename Derived1, typename Derived2>
auto moesp_singular_values(const Eigen::MatrixBase<Derived1>& Y,
                           const Eigen::MatrixBase<Derived2>& U,
                           std::size_t block_rows = 0)
    -> ctrlpp::expected<Eigen::VectorX<typename Derived1::Scalar>, sysid_error>;
```

Returns just the singular values without performing full identification. Useful for model order selection. It refuses the records `moesp` refuses, in the same order, checked at the smallest order `NX = 1`.

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
    if(!sv)
    {
        std::cerr << "MOESP refused the record\n";
        return 1;
    }
    std::cout << "Singular values: " << sv->transpose() << "\n\n";

    // Identify 2nd-order model
    auto result = ctrlpp::moesp<2>(Y, U);
    if(!result)
    {
        std::cerr << "MOESP identification failed\n";
        return 1;
    }

    std::cout << "Identified system:\n"
              << "  A =\n" << result->system.A << "\n"
              << "  B =\n" << result->system.B << "\n"
              << "  C = " << result->system.C << "\n"
              << "  x0 =\n" << result->initial_state << "\n"
              << "  NRMSE = " << result->metrics.nrmse << "\n"
              << "  VAF   = " << result->metrics.vaf << " %\n"
              << "  Cond  = " << result->condition_number << "\n";
}
```

## See Also

- [batch-arx](batch-arx.md)<br/> ARX identification
- [model/state-space](../model/state-space.md)<br/> state-space representation
- [sysid-result](sysid-result.md)<br/> result container types
- [fit-metrics](fit-metrics.md)<br/> NRMSE and VAF metrics
- [guides/sysid/workflow](../../guides/sysid/workflow.md)<br/> system identification workflow guide
- [background/sysid](../../background/sysid.md)<br/> sysid theory and background
