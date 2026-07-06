# particle_filter

Bootstrap Sequential Importance Resampling (SIR) particle filter with compile-time particle count, ESS-adaptive resampling, and roughening to prevent sample degeneracy. Handles arbitrary nonlinear, non-Gaussian state estimation problems by representing the posterior distribution with a weighted set of particles. Supports both log-space and linear weight representations, weighted mean and MAP state extraction, and swappable resampling strategies.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::particle_filter<Scalar, NX, NU, NY, NP, Dynamics, Measurement, Resampler, Rng, Likelihood>` | `#include <ctrlpp/estimation/particle_filter.h>` |

No convenience header exists for this type. Use the categorical path.

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | arithmetic type | Numeric type (e.g. `double`, `float`) |
| `NX` | `std::size_t` | State dimension |
| `NU` | `std::size_t` | Input dimension |
| `NY` | `std::size_t` | Measurement dimension |
| `NP` | `std::size_t` | Number of particles (compile-time) |
| `Dynamics` | satisfies `dynamics_model<Scalar, NX, NU>` | Callable: `(Vector<NX>, Vector<NU>) -> Vector<NX>` |
| `Measurement` | satisfies `measurement_model<Scalar, NX, NY>` | Callable: `(Vector<NX>) -> Vector<NY>` |
| `Resampler` | satisfies `resampling_strategy<Rng, NP>` | Resampling algorithm (default: `systematic_resampling`) |
| `Rng` | satisfies `std::uniform_random_bit_generator` | RNG type (default: `std::mt19937_64`) |
| `Likelihood` | satisfies `pf_likelihood_model<Scalar, NY>` | Measurement-likelihood policy (default: `detail::gaussian_likelihood<Scalar, NY>`) |

## Type Aliases

```cpp
using state_vector_t  = Vector<Scalar, NX>;
using input_vector_t  = Vector<Scalar, NU>;
using output_vector_t = Vector<Scalar, NY>;
```

## Config (`pf_config`)

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `Q` | `Matrix<Scalar, NX, NX>` | Identity | Process noise covariance (used for noise sampling) |
| `R` | `Matrix<Scalar, NY, NY>` | Identity | Measurement noise covariance (used for likelihood) |
| `x0` | `Vector<Scalar, NX>` | Zero | Initial particle mean |
| `P0` | `Matrix<Scalar, NX, NX>` | Identity | Initial particle dispersion covariance |
| `ess_threshold` | `Scalar` | `NP/2` | Effective Sample Size below which resampling triggers |
| `roughening_scale` | `Scalar` | `0.2` | Jitter scale applied after resampling to prevent degeneracy |
| `extraction` | `extraction_method` | `weighted_mean` | State extraction method (`weighted_mean` or `map`) |
| `weights` | `weight_representation` | `log` | Weight storage mode (`log` for numerical stability, `linear` for simplicity) |

## Constructors

```cpp
particle_filter(Dynamics dynamics, Measurement measurement,
                pf_config<Scalar, NX, NU, NY> config, Rng rng = Rng{},
                Likelihood likelihood = Likelihood{});

particle_filter(Dynamics dynamics, Measurement measurement,
                pf_config<Scalar, NX, NU, NY> config, Resampler resampler,
                Rng rng = Rng{}, Likelihood likelihood = Likelihood{});
```

### Factory Function

```cpp
template <std::size_t NP, typename Dynamics, typename Measurement, typename Scalar,
          std::size_t NX, std::size_t NU, std::size_t NY, typename Rng = std::mt19937_64,
          typename Likelihood = detail::gaussian_likelihood<Scalar, NY>>
auto make_particle_filter(Dynamics dynamics, Measurement measurement,
                          pf_config<Scalar, NX, NU, NY> config, Rng rng = Rng{},
                          Likelihood likelihood = Likelihood{});
```

Since NP cannot be deduced via CTAD, this factory function provides a convenient construction interface. Passing a custom `Likelihood` as the trailing argument threads it through to the filter.

## Methods

### predict

```cpp
void predict(const input_vector_t& u);
```

Propagates all NP particles through the dynamics model with additive Gaussian process noise sampled from Q.

### update

```cpp
void update(const output_vector_t& z);
```

Updates particle weights using the configured measurement-likelihood policy (Gaussian by default), normalizes, and resamples (with roughening) if the Effective Sample Size drops below `ess_threshold`. If all particles have negligible likelihood (complete weight collapse), the normalizer resets to uniform weights to recover from particle depletion. This recovery applies to both weight modes: the log mode resets when the maximum log-weight is not finite, and the linear mode resets when the total weight underflows to a non-positive sum.

### state

```cpp
state_vector_t state() const;
```

Returns the state estimate using the configured extraction method (weighted mean or MAP).

### weighted_mean

```cpp
auto weighted_mean() const -> state_vector_t;
```

Returns the weighted mean of all particles.

### map_estimate

```cpp
auto map_estimate() const -> state_vector_t;
```

Returns the particle with the highest weight (maximum a posteriori).

### particles

```cpp
auto particles() const -> const std::array<state_vector_t, NP>&;
```

Returns a const reference to the particle array for inspection or visualization.

## Supporting Types

### systematic_resampling

Low-variance resampling using a single uniform random number and cumulative weight scan. O(NP) complexity.

Header: `#include <ctrlpp/estimation/resampling/systematic_resampling.h>`

### multinomial_resampling

Standard multinomial resampling drawing NP independent samples. Higher variance than systematic.

Header: `#include <ctrlpp/estimation/resampling/multinomial_resampling.h>`

### extraction_method

```cpp
enum class extraction_method { weighted_mean, map };
```

### weight_representation

```cpp
enum class weight_representation { log, linear };
```

Log-space avoids numerical underflow for large particle counts. Linear is simpler but may underflow.

### Likelihood policy

The measurement-likelihood policy is injectable through the last class template parameter and the trailing constructor / factory argument. A policy is a callable that returns a scalar log-likelihood for a single particle's innovation:

```cpp
Scalar operator()(const Vector<Scalar, NY>& z, const Vector<Scalar, NY>& z_pred,
                  const Matrix<Scalar, NY, NY>& R_inv, Scalar log_det_2piR) const;
```

`R_inv` is the precomputed measurement precision and `log_det_2piR` is the precomputed log normalizer `log((2*pi)^NY det R)`, both derived once from the configured `R`. The default policy `detail::gaussian_likelihood<Scalar, NY>` reproduces the Gaussian log-likelihood exactly:

```cpp
return Scalar{-0.5} * (z - z_pred).transpose() * R_inv * (z - z_pred) - Scalar{0.5} * log_det_2piR;
```

A custom policy is the correct place to handle measurement topology. For a bearing sensor, for example, wrap the scalar innovation across the +/-pi cut with `std::remainder(z(0) - z_pred(0), 2*pi)` before forming the Gaussian, so a measurement near +pi and a prediction near -pi are treated as close rather than maximally distant. The policy must satisfy the `pf_likelihood_model<Scalar, NY>` concept.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'true', '' using 1:3 with lines title 'estimate'"

#include <ctrlpp/estimation/particle_filter.h>

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <random>

int main()
{
    using Scalar = double;
    constexpr std::size_t NX = 2, NU = 1, NY = 1, NP = 500;

    auto dynamics = [](const ctrlpp::Vector<Scalar, NX>& x,
                       const ctrlpp::Vector<Scalar, NU>& u) -> ctrlpp::Vector<Scalar, NX> {
        ctrlpp::Vector<Scalar, NX> xn;
        xn(0) = x(0) + 0.1 * x(1);
        xn(1) = x(1) + 0.1 * u(0);
        return xn;
    };

    auto measurement = [](const ctrlpp::Vector<Scalar, NX>& x) -> ctrlpp::Vector<Scalar, NY> {
        ctrlpp::Vector<Scalar, NY> z;
        z(0) = std::atan2(x(1), x(0));  // bearing-only observation
        return z;
    };

    ctrlpp::pf_config<Scalar, NX, NU, NY> cfg{
        .Q = Eigen::Matrix2d::Identity() * 0.01,
        .R = (Eigen::Matrix<Scalar, 1, 1>() << 0.05).finished(),
        .x0 = Eigen::Vector2d::Zero(),
        .P0 = Eigen::Matrix2d::Identity() * 2.0
    };

    auto pf = ctrlpp::make_particle_filter<NP>(dynamics, measurement, cfg, std::mt19937_64{42});

    std::mt19937 rng(123);
    std::normal_distribution<> noise(0.0, std::sqrt(0.05));

    ctrlpp::Vector<Scalar, NX> x_true;
    x_true << 1.0, 0.5;

    for (int k = 0; k < 100; ++k) {
        ctrlpp::Vector<Scalar, NU> u = ctrlpp::Vector<Scalar, NU>::Zero();
        x_true = dynamics(x_true, u);

        ctrlpp::Vector<Scalar, NY> z;
        z(0) = std::atan2(x_true(1), x_true(0)) + noise(rng);

        pf.predict(u);
        pf.update(z);

        auto est = pf.state();
        std::cout << k * 0.1 << "," << x_true(0) << "," << est(0) << "\n";
    }
}
```

## See Also

- [ukf](ukf.md)<br/> sigma-point filter for moderate nonlinearity
- [ekf](ekf.md)<br/> linearization-based filter for smooth systems
- [observer-policy](observer-policy.md)<br/> concept satisfied by this type
- [background/particle-filter](../../background/particle-filter.md)<br/> SIR derivation
- [guides/estimation/observer-controller](../../guides/estimation/observer-controller.md)<br/> composing observers with controllers
