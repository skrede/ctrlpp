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

## Construction

```cpp
static auto create(Dynamics dynamics, Measurement measurement,
                   pf_config<Scalar, NX, NU, NY> config, Rng rng = Rng{},
                   Likelihood likelihood = Likelihood{})
    -> ctrlpp::expected<particle_filter, filter_error>;

static auto create(Dynamics dynamics, Measurement measurement,
                   pf_config<Scalar, NX, NU, NY> config, Resampler resampler,
                   Rng rng = Rng{}, Likelihood likelihood = Likelihood{})
    -> ctrlpp::expected<particle_filter, filter_error>;
```

`create` is the only construction path. The second overload takes the
resampling strategy; the first default-constructs it. The constructor is public
only so that the result can be built in place, and it takes a passkey that only
the filter can make, so no caller can reach it without the factory.

Every refusal happens before any filter state exists. On success the filter is
constructed directly in the storage the caller initializes from the result, and
guaranteed copy elision carries that through `make_particle_filter` as well, so
the factory's frame never holds a whole filter. A `static` or a class member
initialized from `create` keeps the filter off a task stack altogether:

```cpp
static auto filter = make_particle_filter<4000>(dynamics, measurement, config);
if(!filter)
    return;
filter->predict(u);
```

Measured with `g++ (GCC) 16.2.1`, `-std=c++20 -O2 -fno-exceptions -fno-rtti
-fstack-usage`, Eigen 3.4.1, x86-64 Linux, `double`, NX 2, NY 1, the default
resampler, generator and likelihood, the factory called through a non-inlined
wrapper into a function-local `static`:

| NP | `sizeof(particle_filter)` | Factory frames | `resample_log` frame | `resample_linear` frame |
|---:|---:|---:|---:|---:|
| 100 | 8,256 | 8,352 | 2,464 | 1,648 |
| 1,000 | 58,656 | 8,352 | 24,064 | 16,048 |
| 4,000 | 226,656 | 8,352 | 96,064 | 64,048 |

The factory's frames do not depend on NP. The object holds the particles twice
(the set and the set before the current resampling) and the weights three times
(log, linear and staged). The resampling's frame does grow with NP: its index
array, the resampler's cumulative sum and, in log mode, the linear weights are
locals of NP entries each.

```cpp
template <std::size_t NP, typename Dynamics, typename Measurement, typename Scalar,
          std::size_t NX, std::size_t NU, std::size_t NY, typename Rng = std::mt19937_64,
          typename Likelihood = detail::gaussian_likelihood<Scalar, NY>>
auto make_particle_filter(Dynamics dynamics, Measurement measurement,
                          pf_config<Scalar, NX, NU, NY> config, Rng rng = Rng{},
                          Likelihood likelihood = Likelihood{})
    -> ctrlpp::expected<particle_filter<Scalar, NX, NU, NY, NP, Dynamics, Measurement,
                                        systematic_resampling, Rng, Likelihood>,
                        filter_error>;
```

`make_particle_filter` is `create` with the particle count named first, since NP
cannot be deduced from the arguments. Passing a custom `Likelihood` as the
trailing argument threads it through to the filter.

Rejections, checked in the order `pf_config` declares its fields, and for each
matrix in the order finiteness, symmetry, definiteness (then, for R,
invertibility). The bound for an N x N matrix M is
`max(N eps max|M_ij|, N denorm_min)`. As a definiteness floor it is the order of
the symmetric eigensolver's backward error, below which a computed eigenvalue is
indistinguishable from zero, kept from underflowing to zero for a
subnormal-scale matrix. An asymmetry `max|M_ij - M_ji|` within it is rounding
drift: the matrix is not refused but symmetrized (each differing pair replaced
by its midpoint), and the definiteness verdict, the sampling factor, R's
inverse and its log-determinant are all formed from that symmetrized matrix.
The convention and its count are in
[Rounding drift is accepted, structure is refused](../../guides/patterns/numerical-behavior.md#rounding-drift-is-accepted-structure-is-refused).

| Enumerator | Condition | Why it is refused |
|------------|-----------|-------------------|
| `non_finite_process_noise` | Q has a non-finite entry | No Gaussian to sample the process noise from |
| `asymmetric_process_noise` | Q's asymmetry exceeds the bound | Not a covariance, so no Gaussian to sample |
| `indefinite_process_noise` | Q has an eigenvalue below `-floor` | Not a covariance, so no Gaussian to sample |
| `non_finite_measurement_noise` | R has a non-finite entry | No density to score a particle with |
| `asymmetric_measurement_noise` | R's asymmetry exceeds the bound | Not a covariance, so no density to score with |
| `indefinite_measurement_noise` | R has an eigenvalue below `-floor` | Not a covariance, so no density to score with |
| `singular_measurement_noise` | `lambda_min(R) <= NY eps lambda_max(R)` | R is inverted and its log-determinant enters every weight; a singular R has neither |
| `non_finite_initial_state` | x0 has a non-finite component | Every particle is drawn around it |
| `non_finite_initial_covariance` | P0 has a non-finite entry | No Gaussian prior to sample |
| `asymmetric_initial_covariance` | P0's asymmetry exceeds the bound | Not a covariance, so no Gaussian prior |
| `indefinite_initial_covariance` | P0 has an eigenvalue below `-floor` | Not a covariance, so no Gaussian prior |

Definiteness is a domain condition for this filter, not a numerical-behavior
question: it samples from N(0, Q) and N(x0, P0) and scores with R's inverse and
log-determinant. A singular but positive semidefinite Q or P0 is accepted, since
a rank-one process noise (the discretized white-noise-acceleration model) and a
point-mass prior are legitimate. The sampling factor is the Cholesky factor when
the matrix is positive definite, and otherwise `V diag(sqrt(max(lambda, 0)))`
from the symmetric eigendecomposition, which keeps every zero eigenvalue at zero
so the particles carry exactly the covariance they were given.

## Methods

### predict

```cpp
void predict(const input_vector_t& u);
```

Propagates all NP particles through the dynamics model with additive Gaussian process noise sampled from Q.

### update

```cpp
auto update(const output_vector_t& z)
    -> ctrlpp::expected<void, pf_update_error>;
```

The step is refused **before any weight is committed and before the random generator is drawn from**, so a refused step leaves the particles, the weights and the generator bitwise unchanged and the caller may continue with the next sample. The checks run in this order, the carried estimate ahead of the measurement:

| Condition | Error | Why it is refused |
|---|---|---|
| a carried particle is already non-finite | `pf_update_error::non_finite_state` | No weight computed from it means anything; `health()` latches `non_finite_estimate` |
| the supplied measurement has a non-finite component | `pf_update_error::non_finite_measurement` | Every weight would be non-finite and the posterior would be discarded |
| the measurement is finite but a new weight would be NaN or +inf | `pf_update_error::non_finite_result` | A custom likelihood policy or a measurement model produced a non-finite likelihood |

The new weights are computed into a staging array and committed only when none is NaN or +inf; a log-weight of -inf (a linear weight of zero) is a legitimate zero likelihood. On a step that runs, the weights are updated with the configured measurement-likelihood policy (Gaussian by default), normalized, and the particles are resampled (with roughening) if the Effective Sample Size drops below `ess_threshold`.

When no particle explains a finite measurement (every likelihood is zero), the normalizer resets to uniform weights and `health()` reports `pf_health::weights_reset`. The step still succeeds: the posterior is valid but uninformed. This applies to both weight modes: the log mode resets when the maximum log-weight is -inf, and the linear mode resets when the total weight underflows to zero.

After every resampling the particles are checked for finiteness. Roughening scales its jitter with the cloud's coordinate range, so a cloud near the edge of the representable range can overflow it: the range itself overflows to infinity, or the jitter carries a particle past the largest finite value. The resampling is then undone: the particles are restored to the set before it, the step's conditioned weights are kept (they are the posterior over that set, whose particles passed the step's own finiteness check), the step succeeds, and `health()` latches `pf_health::resampling_rolled_back`. So `update` never reports success with a non-finite particle. The generator has advanced by the draws the undone resampling made.

`predict` is not fallible. A prediction that makes the particles non-finite is reported by the next `update` as `non_finite_state` and through [`health`](#health).

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

### covariance

```cpp
auto covariance() const -> Matrix<Scalar, NX, NX>;
```

Returns the weighted second central moment of the particle set about its weighted mean: the posterior covariance the cloud and its weights represent. This is the uncertainty half of the estimate, and it is not recoverable from `particles()` alone, because the unweighted dispersion of the array ignores the weights and therefore describes the prior spread whenever the last `update` did not trigger a resampling.

Three properties of the definition:

- The center is always the weighted mean, including when `extraction` is `map`. A second moment about any other point is larger than the covariance and is not one; dispersion about the MAP estimate can be formed from `particles()` and `map_estimate()`.
- There is no Bessel correction. The weights sum to one, so this is the weighted second moment, matching the convention the unscented filter's sigma-point covariance uses.
- Uniform weights need no special case. The expression reduces exactly to the plain second moment of the particles about their plain mean, which is the correct answer after a depletion that `health()` reports as `weights_reset`: the measurement carried no information, so the reported uncertainty is the dispersion the filter was already carrying.

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

### health

```cpp
auto health() const -> pf_health;
```

Returns the persistent health status. The enumerators are ordered by severity and the status latches, never returning to a lower rung on its own.

| Enumerator | Meaning |
|---|---|
| `pf_health::ok` | Every update so far conditioned the particles on its measurement |
| `pf_health::weights_reset` | A finite measurement that no particle explains reset the weights to uniform; the step succeeded with an uninformed posterior |
| `pf_health::resampling_rolled_back` | A resampling whose roughened particles would have been non-finite was undone; the step succeeded with its weighted posterior over the set before the resampling |
| `pf_health::non_finite_estimate` | A step found the particles already non-finite, which is how a poisoned control or dynamics model becomes visible |

A refused measurement does not set it: the refusal mutates nothing.

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

    auto created = ctrlpp::make_particle_filter<NP>(dynamics, measurement, cfg, std::mt19937_64{42});
    if (!created.has_value()) {
        std::cerr << "invalid particle filter configuration\n";
        return 1;
    }
    auto& pf = *created;

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
        if (const auto stepped = pf.update(z); !stepped.has_value()) {
            std::cerr << "measurement refused at step " << k << "\n";
            continue;
        }

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
