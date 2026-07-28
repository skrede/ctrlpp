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
| `lambda` | `Scalar` | `0.99` | Forgetting factor, `0 < lambda <= 1`. Lower values forget faster. |
| `P0` | `Matrix<Scalar, NP, NP>` | `1000 * I` | Initial covariance matrix |
| `cov_upper_bound` | `Scalar` | `1e6` | Maximum trace-per-dimension before covariance is scaled down |

## Construction

```cpp
static auto create(rls_config<Scalar, NP> config = {})
    -> ctrlpp::expected<rls, rls_error>;
```

`create` is the only construction path; there is no public constructor. A
factory that sat beside one would validate nothing, since any caller could
bypass it. Parameters are initialized to zero and the covariance to `P0`.

| Enumerator | Condition | Why |
|------------|-----------|-----|
| `rls_error::non_positive_forgetting_factor` | `lambda` is not finite and strictly positive | Read off the covariance update, where `lambda` is the **divisor**: `P <- (P - k phi^T P) / lambda`. At zero the update divides by zero and `P` is non-finite from the first sample; at a negative value the division negates a positive semidefinite matrix, so every later gain points against the error rather than along it. This is an exact domain condition, not a preference. |
| `rls_error::forgetting_factor_above_unity` | `lambda > 1` | A **contract** violation, not a domain violation, and it carries its own enumerator so the two are not conflated. The arithmetic here is well defined: dividing by a factor above one deflates the covariance faster than the measurement update alone, so the gain collapses toward zero and the estimator silently stops adapting -- it does not diverge. Exponential forgetting is defined on `(0, 1]`, which is what this rejection enforces. `lambda == 1` is accepted: that is ordinary recursive least squares with no forgetting. |
| `rls_error::non_finite_initial_covariance` | `P0` has a non-finite entry | `P0` seeds the recursion, which has no mechanism that returns a non-finite covariance to a finite one. Finiteness is the whole test: a **singular or zero `P0` is accepted**, being a legitimate starting point that says only "I am certain of these parameters". |
| `rls_error::non_positive_covariance_bound` | `cov_upper_bound` is not finite and strictly positive | Once the trace exceeds `NP * cov_upper_bound` the update rescales `P` by `bound*NP/trace`, so a negative bound negates the covariance and a zero bound drives it to exactly zero on the first sample, after which the gain is zero forever and the parameters never move again. A non-finite bound disables the clamp the field exists to impose. |

```cpp
auto estimator = ctrlpp::rls<double, NP>::create(cfg);
if(!estimator)
    return estimator.error();
```

## Methods

### update

```cpp
ctrlpp::expected<void, rls_update_error> update(Scalar y, const Vector<Scalar, NP>& phi);
```

Incorporates a new observation. Given measurement `y` and regressor vector `phi`, updates the parameter estimate and covariance using the standard RLS gain computation with forgetting factor.

The cycle is classified before any member is written, so a refused cycle leaves the parameters and the covariance **bitwise unchanged** and the caller may retry on the next sample. That matters here more than a boolean suggested: those two members are the estimator's entire memory, nothing re-derives them, so one admitted non-finite sample makes them non-finite forever while `parameters()` keeps returning a vector the caller has no way to distrust.

| Enumerator | Condition | Why it has its own name |
|---|---|---|
| `rls_update_error::non_finite_state` | the carried covariance or parameter vector is already non-finite | Nothing this cycle can produce is meaningful, and the recursion has no mechanism that returns a non-finite covariance to a finite one, so only reconstruction recovers. |
| `rls_update_error::non_finite_observation` | `y` is not finite | It enters the prediction error and from there the parameter step directly. The repair is upstream, in whatever measures the output. |
| `rls_update_error::non_finite_regressor` | `phi` has a non-finite component | It reaches the parameters through the gain and the covariance through the rank-one update. A different subsystem from the one that measures the output, hence a different enumerator. |
| `rls_update_error::non_finite_denominator` | every operand was finite and `lambda + phi^T P phi` still left the scalar's range | A magnitude fault rather than a domain violation of either operand; the repair is a rescaling, not a replacement. |
| `rls_update_error::indefinite_covariance` | the denominator is negative by more than its own resolution | `lambda + phi^T P phi` with `lambda > 0` cannot be negative for a positive semidefinite `P`, so this says `P` stopped being a covariance. A gain formed from it points **against** the prediction error. |
| `rls_update_error::denominator_below_resolution` | the denominator carries no significant digits at the scale of the operands that formed it | Two situations reach it and it deliberately does **not** distinguish them, because the arithmetic cannot: a regressor the covariance genuinely cannot resolve, and a cancellation in `phi^T P phi` leaving a residue made entirely of rounding. |

The resolution floor is **derived, not chosen**: it is the counted rounding of the two contractions that formed the denominator -- `2 * (2 * NP - 1) + 1` operations -- times the scalar's machine epsilon, times the largest operand that entered, which is the forgetting factor or the Cauchy-Schwarz bound on the quadratic form. An absolute floor is wrong in both directions and both are reachable: on a problem posed far below unit scale it refuses a perfectly well conditioned update, and on one posed far above it accepts a denominator whose significant digits have all cancelled away.

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
        if(const auto applied = estimator.update(y, phi); !applied)
        {
            std::cerr << "RLS refused sample " << k << "\n";
            return 1;
        }

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

- [recursive-arx](recursive-arx.md)<br/> recursive ARX identification using RLS
- [batch-arx](batch-arx.md)<br/> batch ARX identification
- [guides/sysid/workflow](../../guides/sysid/workflow.md)<br/> system identification workflow guide
- [background/sysid](../../background/sysid.md)<br/> sysid theory and background
