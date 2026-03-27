# rls

Recursive Least Squares estimator with exponential forgetting factor and bounded covariance. Suitable for online parameter estimation where measurements arrive one at a time. The forgetting factor controls the effective memory length, allowing the estimator to track slowly time-varying parameters.

## Header and Alias

| Form | Header |
|------|--------|
| `rls<Scalar, NP>` | `#include <ctrlpp/sysid/rls.h>` |
| (convenience) | `#include <ctrlpp/sysid.h>` |

```cpp
template <typename Scalar, std::size_t NP>
class rls;
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | floating-point | Numeric type (`double`, `float`) |
| `NP` | `>= 1` | Number of parameters to estimate |

## rls_config

Configuration struct `rls_config<Scalar, NP>` passed at construction.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `lambda` | `Scalar` | `0.99` | Forgetting factor (0 < lambda <= 1). Lower values forget faster. |
| `P0` | `Matrix<Scalar, NP, NP>` | `1000 * I` | Initial covariance matrix |
| `cov_upper_bound` | `Scalar` | `1e6` | Maximum trace-per-dimension before covariance is scaled down |

## Constructors

```cpp
explicit rls(rls_config<Scalar, NP> config = {});
```

Constructs the estimator from configuration. Parameters initialised to zero, covariance to `P0`.

## Methods

### update

```cpp
void update(Scalar y, const Vector<Scalar, NP>& phi);
```

Incorporates a new observation. Given measurement `y` and regressor vector `phi`, updates the parameter estimate and covariance using the standard RLS gain computation with forgetting factor.

### parameters

```cpp
const Vector<Scalar, NP>& parameters() const;
```

Returns the current parameter estimate vector.

### covariance

```cpp
const Matrix<Scalar, NP, NP>& covariance() const;
```

Returns the current covariance matrix.

## Usage Example

```cpp
#include <ctrlpp/sysid/rls.h>

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <random>

int main()
{
    // Identify a first-order system: y(t) = a*y(t-1) + b*u(t-1)
    constexpr std::size_t NP = 2;

    ctrlpp::rls_config<double, NP> cfg{
        .lambda = 0.98,
        .P0 = Eigen::Matrix2d::Identity() * 100.0};

    ctrlpp::rls<double, NP> estimator(cfg);

    // True system: y(t) = 0.8*y(t-1) + 0.5*u(t-1) + noise
    constexpr double a_true = 0.8;
    constexpr double b_true = 0.5;

    std::mt19937 rng(42);
    std::normal_distribution<double> noise(0.0, 0.01);
    std::normal_distribution<double> input_dist(0.0, 1.0);

    double y_prev = 0.0;
    double u_prev = 0.0;

    for(int k = 0; k < 200; ++k)
    {
        double u = input_dist(rng);
        double y = a_true * y_prev + b_true * u_prev + noise(rng);

        Eigen::Vector2d phi(y_prev, u_prev);
        estimator.update(y, phi);

        if(k % 50 == 49)
        {
            auto theta = estimator.parameters();
            std::cout << "k=" << k
                      << "  a_hat=" << theta[0] << " (true=" << a_true << ")"
                      << "  b_hat=" << theta[1] << " (true=" << b_true << ")\n";
        }

        y_prev = y;
        u_prev = u;
    }
}
```

## See Also

- [recursive-arx](recursive-arx.md) -- recursive ARX identification using RLS
- [batch-arx](batch-arx.md) -- batch ARX identification
- [guides/sysid/workflow](../guides/sysid/workflow.md) -- system identification workflow guide
- [reference/sysid-theory](../reference/sysid-theory.md) -- sysid theory and background
