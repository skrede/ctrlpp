# recursive_arx

Recursive ARX model identification using RLS internally. Processes input/output data one sample at a time and maintains a running ARX(NA, NB) parameter estimate. Can convert the current estimate to an observer canonical form discrete state-space model at any time.

## Header and Alias

| Form | Header |
|------|--------|
| `recursive_arx<Scalar, NA, NB, NU, NY>` | `#include <ctrlpp/sysid/recursive_arx.h>` |
| (convenience) | `#include <ctrlpp/sysid.h>` |

```cpp
template <typename Scalar, std::size_t NA, std::size_t NB,
          std::size_t NU = 1, std::size_t NY = 1>
class recursive_arx;
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | floating-point | Numeric type |
| `NA` | `>= 1` | Number of auto-regressive (output) terms |
| `NB` | `>= 1` | Number of exogenous (input) terms |
| `NU` | `>= 1` | Number of inputs (default 1) |
| `NY` | `>= 1` | Number of outputs (default 1) |

## Constructors

```cpp
explicit recursive_arx(rls_config<Scalar, NP> config = {});
```

Where `NP = NA * NY + NB * NU`. Constructs the identifier with optional RLS configuration (forgetting factor, initial covariance, covariance bound).

## Methods

### update

```cpp
void update(Scalar y, Scalar u);
```

Processes a new input/output pair. Builds the regressor vector from the internal history buffers and updates the RLS parameter estimate.

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
discrete_state_space<Scalar, NA, NU, NY> to_state_space() const;
```

Converts the current parameter estimate to observer canonical form state-space matrices (A, B, C, D).

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

    ctrlpp::recursive_arx<double, NA, NB> identifier;

    std::mt19937 rng(42);
    std::normal_distribution<double> input_dist(0.0, 1.0);
    std::normal_distribution<double> noise(0.0, 0.01);

    double y_prev = 0.0;
    double y_prev2 = 0.0;

    for(int k = 0; k < 300; ++k)
    {
        double u = input_dist(rng);
        double y = 0.7 * y_prev - 0.2 * y_prev2 + 0.5 * u + noise(rng);

        identifier.update(y, u);

        if(k % 100 == 99)
        {
            auto theta = identifier.parameters();
            std::cout << "k=" << k << "  params=[" << theta.transpose() << "]\n";

            auto sys = identifier.to_state_space();
            std::cout << "  A=\n" << sys.A << "\n  B=\n" << sys.B << "\n";
        }

        y_prev2 = y_prev;
        y_prev = y;
    }
}
```

## See Also

- [batch-arx](batch-arx.md) -- offline batch ARX identification
- [rls](rls.md) -- underlying recursive least squares estimator
