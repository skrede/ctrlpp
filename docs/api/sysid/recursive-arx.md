# recursive_arx

Recursive ARX model identification using RLS internally. Processes input/output data one sample at a time and maintains a running ARX(NA, NB) parameter estimate. Can convert the current estimate to an observer canonical form discrete state-space model at any time.

The estimator is single-input single-output: `update` takes one scalar output and one scalar input, and the realization it produces has one input and one output. Multi-channel (MIMO) ARX identification is not implemented, and there is no way to ask this type for it.

## Header and Alias

| Form | Header |
|------|--------|
| `recursive_arx<Scalar, NA, NB>` | `#include <ctrlpp/sysid/recursive_arx.h>` |
| (convenience) | `#include <ctrlpp/sysid.h>` |

```cpp
template <typename Scalar, std::size_t NA, std::size_t NB>
class recursive_arx;
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | floating-point | Numeric type |
| `NA` | `>= 1` | Number of auto-regressive (output) terms |
| `NB` | `>= 1` | Number of exogenous (input) terms |

## Construction

```cpp
static auto create(rls_config<Scalar, NP> config = {})
    -> ctrlpp::expected<recursive_arx, rls_error>;
```

Where `NP = NA + NB`. `create` is the only construction path; there is
no public constructor. It builds the identifier from an optional RLS
configuration (forgetting factor, initial covariance, covariance bound).

This type owns a `rls` instance and configures it from the same aggregate, so it
has no configuration condition of its own: it **forwards** that estimator's
rejection verbatim rather than restating the conditions, where a second copy
could drift out of step with the arithmetic it describes. See
[`rls`](rls.md#construction) for the enumerators and their derivations.

```cpp
auto identifier = ctrlpp::recursive_arx<double, NA, NB>::create(cfg);
if(!identifier)
    return identifier.error();
```

## Methods

### update

```cpp
ctrlpp::expected<void, recursive_arx_update_error> update(Scalar y, Scalar u);
```

Processes a new input/output pair. Builds the regressor vector from the internal history buffers and updates the RLS parameter estimate.

The underlying estimator's refusal is **forwarded verbatim**: every enumerator of [`rls_update_error`](rls.md#update) has a same-named counterpart in `recursive_arx_update_error`, so the reason is never restated under a name that could drift out of step with the arithmetic it describes. One enumerator is the wrapper's own: `non_finite_input`, reported when `u` is not finite, before the regressor is built.

A refused cycle also leaves the regressor history, the write index and the sample count untouched -- recording a sample the estimator refused as non-finite would poison every regressor built afterwards, so the poison would latch here after the estimator had correctly declined it.

### parameters

```cpp
const Vector<Scalar, NP>& parameters() const;
```

Returns the current raw parameter vector `[a1, ..., aNa, b1, ..., bNb]`.

### covariance

```cpp
const Matrix<Scalar, NP, NP>& covariance() const;
```

Returns the current RLS covariance matrix.

### to_state_space

```cpp
discrete_state_space<Scalar, std::max(NA, NB), 1, 1> to_state_space() const;
```

Converts the current parameter estimate to observer canonical form state-space matrices (A, B, C, D). The realization has `max(NA, NB)` states so that every b-coefficient is represented even when `NB > NA` (Ljung 1999, Ch. 4), one input and one output.

## Usage Example

```cpp
#include <ctrlpp/sysid/recursive_arx.h>

#include <Eigen/Dense>

#include <iostream>
#include <random>

int main()
{
    // Online identification of y(t) = 0.7*y(t-1) - 0.2*y(t-2) + 0.5*u(t-1)
    constexpr std::size_t NA = 2;
    constexpr std::size_t NB = 1;

    auto created = ctrlpp::recursive_arx<double, NA, NB>::create();
    if(!created)
    {
        std::cerr << "recursive ARX refused its configuration\n";
        return 1;
    }
    auto& identifier = *created;

    std::mt19937 rng(42);
    std::normal_distribution<double> input_dist(0.0, 1.0);
    std::normal_distribution<double> noise(0.0, 0.01);

    double y_prev = 0.0;
    double y_prev2 = 0.0;
    double u_prev = 0.0;

    for(int k = 0; k < 300; ++k)
    {
        double u = input_dist(rng);
        double y = 0.7 * y_prev - 0.2 * y_prev2 + 0.5 * u_prev + noise(rng);

        if(const auto applied = identifier.update(y, u); !applied)
        {
            std::cerr << "recursive ARX refused sample " << k << "\n";
            return 1;
        }

        if(k % 100 == 99)
        {
            auto theta = identifier.parameters();
            std::cout << "k=" << k << "  params=[" << theta.transpose() << "]\n";

            auto sys = identifier.to_state_space();
            std::cout << "  A=\n" << sys.A << "\n  B=\n" << sys.B << "\n";
        }

        u_prev = u;
        y_prev2 = y_prev;
        y_prev = y;
    }
}
```

## See Also

- [batch-arx](batch-arx.md)<br/> offline batch ARX identification
- [rls](rls.md)<br/> underlying recursive least squares estimator
- [sysid-result](sysid-result.md)<br/> result container types
- [guides/sysid/workflow](../../guides/sysid/workflow.md)<br/> system identification workflow guide
- [background/sysid](../../background/sysid.md)<br/> sysid theory and background
