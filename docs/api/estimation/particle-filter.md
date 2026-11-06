# particle_filter

Bootstrap Sequential Importance Resampling (SIR) particle filter with compile-time particle count, ESS-adaptive resampling, and roughening to prevent sample degeneracy. It represents the posterior with a weighted set of particles, so the posterior need not be Gaussian and the models may be strongly nonlinear. The process noise is additive Gaussian with covariance `Q`, and the measurement likelihood is Gaussian in `R` unless a likelihood policy is injected (see [Limitations](#limitations)). Supports both log-space and linear weight representations, weighted mean and highest-weight state extraction, and swappable resampling strategies.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::particle_filter<Scalar, NX, NU, NY, NP, Dynamics, Measurement, Resampler, Rng, Likelihood>` | `#include <ctrlpp/estimation/particle_filter.h>` |
| `ctrlpp::particle_filter<...>`, `ctrlpp::make_particle_filter` | `#include <ctrlpp/particle_filter.h>` (convenience) |
| `ctrlpp::pf_config`, `ctrlpp::extraction_method`, `ctrlpp::weight_representation`, `ctrlpp::pf_update_error`, `ctrlpp::pf_health`, `ctrlpp::pf_likelihood_model` | `#include <ctrlpp/estimation/particle_filter_types.h>` (included by `particle_filter.h`) |
| `ctrlpp::resampling_strategy` | `#include <ctrlpp/estimation/resampling/resampling_strategy.h>` (included by `particle_filter.h`) |
| `ctrlpp::systematic_resampling` | `#include <ctrlpp/estimation/resampling/systematic_resampling.h>` (included by `particle_filter.h`) |
| `ctrlpp::multinomial_resampling` | `#include <ctrlpp/estimation/resampling/multinomial_resampling.h>` (not included by `particle_filter.h`) |

`particle_filter` satisfies `ObserverPolicy` ([observer-policy](observer-policy.md)). It is not a `CovarianceObserver`: it has `covariance()` but no `innovation()`.

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | `ctrlpp_floating_scalar` | Floating-point type (e.g. `double`, `float`) |
| `NX` | `std::size_t`, `NX > 0` | State dimension |
| `NU` | `std::size_t`, `NU > 0` | Input dimension |
| `NY` | `std::size_t`, `NY > 0` | Measurement dimension |
| `NP` | `std::size_t`, `NP > 0` | Number of particles (compile-time) |
| `Dynamics` | satisfies `dynamics_model<Dynamics, Scalar, NX, NU>` | Callable: `(Vector<Scalar, NX>, Vector<Scalar, NU>) -> Vector<Scalar, NX>` |
| `Measurement` | satisfies `measurement_model<Measurement, Scalar, NX, NY>` | Callable: `(Vector<Scalar, NX>) -> Vector<Scalar, NY>` |
| `Resampler` | satisfies `resampling_strategy<Resampler, Rng, NP>` | Resampling algorithm (default: `systematic_resampling`) |
| `Rng` | satisfies `std::uniform_random_bit_generator` | RNG type (default: `std::mt19937_64`) |
| `Likelihood` | satisfies `pf_likelihood_model<Likelihood, Scalar, NY>` | Measurement-likelihood policy (default: `detail::gaussian_likelihood<Scalar, NY>`) |

## Type Aliases

```cpp
using state_vector_t  = Vector<Scalar, NX>;
using input_vector_t  = Vector<Scalar, NU>;
using output_vector_t = Vector<Scalar, NY>;
```

## Config (`pf_config`)

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `Q` | `Matrix<Scalar, NX, NX>` | Identity | Discrete process noise covariance, sampled once per particle per `predict` call |
| `R` | `Matrix<Scalar, NY, NY>` | Identity | Measurement noise covariance, handed to the likelihood policy as its inverse and `log((2 pi)^NY det R)` |
| `x0` | `Vector<Scalar, NX>` | Zero | Initial particle mean |
| `P0` | `Matrix<Scalar, NX, NX>` | Identity | Initial particle dispersion covariance |
| `ess_threshold` | `Scalar` | `-1` | Effective Sample Size below which resampling triggers; a negative value selects `NP / 2` |
| `roughening_scale` | `Scalar` | `0.2` | Roughening constant applied after resampling; zero or negative disables roughening |
| `extraction` | `extraction_method` | `weighted_mean` | What `state()` returns (`weighted_mean` or `map`) |
| `weights` | `weight_representation` | `log` | Weight storage mode (`log` for numerical stability, `linear` for simplicity) |

## Construction

```cpp
static ctrlpp::expected<particle_filter, filter_error> create(
    Dynamics dynamics, Measurement measurement, pf_config<Scalar, NX, NU, NY> config,
    Rng rng = Rng{}, Likelihood likelihood = Likelihood{});

static ctrlpp::expected<particle_filter, filter_error> create(
    Dynamics dynamics, Measurement measurement, pf_config<Scalar, NX, NU, NY> config,
    Resampler resampler, Rng rng = Rng{}, Likelihood likelihood = Likelihood{});
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
static auto filter = ctrlpp::make_particle_filter<4000>(dynamics, measurement, config);
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
template <std::size_t NP, typename Dynamics, typename Measurement, ctrlpp_floating_scalar Scalar,
          std::size_t NX, std::size_t NU, std::size_t NY, typename Rng = std::mt19937_64,
          typename Likelihood = detail::gaussian_likelihood<Scalar, NY>>
ctrlpp::expected<particle_filter<Scalar, NX, NU, NY, NP, Dynamics, Measurement,
                                 systematic_resampling, Rng, Likelihood>,
                 filter_error>
make_particle_filter(Dynamics dynamics, Measurement measurement,
                     pf_config<Scalar, NX, NU, NY> config, Rng rng = Rng{},
                     Likelihood likelihood = Likelihood{});
```

`make_particle_filter` is `create` with the particle count named first, since NP
cannot be deduced from the arguments. It always uses `systematic_resampling`; a
different resampler goes through the second `create` overload. Passing a custom
`Likelihood` as the trailing argument threads it through to the filter.

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
so the particles carry exactly the covariance they were given. `ess_threshold`,
`roughening_scale`, `extraction` and `weights` are not validated. The particles
are drawn from N(x0, P0) at construction, with uniform weights.

## Methods

### predict

```cpp
void predict(const input_vector_t& u);
```

Propagates every particle through the dynamics model and adds an independent draw from N(0, Q): `x_i = f(x_i, u) + w_i`. `Q` is added once per call, unscaled by any time step. The weights are unchanged.

### update

```cpp
ctrlpp::expected<void, pf_update_error> update(const output_vector_t& z);
```

The step is refused **before any weight is committed and before the random generator is drawn from**, so a refused step leaves the particles, the weights and the generator bitwise unchanged and the caller may continue with the next sample. The checks run in this order, the carried estimate ahead of the measurement:

| Condition | Error | Why it is refused |
|---|---|---|
| a carried particle is already non-finite | `pf_update_error::non_finite_state` | No weight computed from it means anything; `health()` latches `non_finite_estimate` |
| the supplied measurement has a non-finite component | `pf_update_error::non_finite_measurement` | Every weight would be non-finite and the posterior would be discarded |
| the measurement is finite but a new weight would be NaN or +inf | `pf_update_error::non_finite_result` | A custom likelihood policy or a measurement model produced a non-finite likelihood |

The new weights are computed into a staging array and committed only when none is NaN or +inf; a log-weight of -inf (a linear weight of zero) is a legitimate zero likelihood. On a step that runs, each weight is multiplied by the particle's likelihood from the configured policy (Gaussian by default; in log mode the log-likelihood is added), the weights are normalized (log-sum-exp in log mode), and the particles are resampled if the Effective Sample Size `1 / sum_i w_i^2` is below `ess_threshold`. Without a resampling the weights carry into the next update. A resampling draws NP indices with the resampler, replaces the set by the indexed copies, roughens it, and resets the weights to `1 / NP`.

Roughening jitters each coordinate `d` independently with N(0, sigma_d^2), `sigma_d = roughening_scale * E_d * NP^(-1/NX)`, where `E_d` is the range (maximum minus minimum) of that coordinate over the resampled set. A coordinate whose range is zero falls back to `sigma_d = roughening_scale * sqrt(Q_dd)`, and a coordinate whose sigma is still zero is left unjittered.

When no particle explains a finite measurement (every likelihood is zero), the normalizer resets to uniform weights and `health()` reports `pf_health::weights_reset`. The step still succeeds: the posterior is valid but uninformed. This applies to both weight modes: the log mode resets when the maximum log-weight is -inf, and the linear mode resets when the total weight underflows to zero.

After every resampling the particles are checked for finiteness. Roughening scales its jitter with the cloud's coordinate range, so a cloud near the edge of the representable range can overflow it: the range itself overflows to infinity, or the jitter carries a particle past the largest finite value. The resampling is then undone: the particles are restored to the set before it, the step's conditioned weights are kept (they are the posterior over that set, whose particles passed the step's own finiteness check), the step succeeds, and `health()` latches `pf_health::resampling_rolled_back`. So `update` never reports success with a non-finite particle. The generator has advanced by the draws the undone resampling made.

`predict` is not fallible. A prediction that makes the particles non-finite is reported by the next `update` as `non_finite_state` and through [`health`](#health). A refusal with `non_finite_state` latches `non_finite_estimate`; one with `non_finite_measurement` or `non_finite_result` leaves `health()` as it was.

### state

```cpp
state_vector_t state() const;
```

Returns `weighted_mean()`, or `map_estimate()` when `extraction` is `extraction_method::map`.

### weighted_mean

```cpp
state_vector_t weighted_mean() const;
```

Returns the weighted mean of all particles.

### covariance

```cpp
Matrix<Scalar, NX, NX> covariance() const;
```

Returns the weighted second central moment of the particle set about its weighted mean: the posterior covariance the cloud and its weights represent. This is the uncertainty half of the estimate, and it is not recoverable from `particles()` alone, because the unweighted dispersion of the array ignores the weights and therefore describes the prior spread whenever the last `update` did not trigger a resampling.

Three properties of the definition:

- The center is always the weighted mean, including when `extraction` is `map`. A second moment about any other point is larger than the covariance and is not one; dispersion about the MAP estimate can be formed from `particles()` and `map_estimate()`.
- There is no Bessel correction. The weights sum to one, so this is the weighted second moment, matching the convention the unscented filter's sigma-point covariance uses.
- Uniform weights need no special case. The expression reduces exactly to the plain second moment of the particles about their plain mean, which is the correct answer after a depletion that `health()` reports as `weights_reset`: the measurement carried no information, so the reported uncertainty is the dispersion the filter was already carrying.

### map_estimate

```cpp
state_vector_t map_estimate() const;
```

Returns the particle with the highest weight, the lowest index on a tie. This is the highest-weight sample, not an estimate of the posterior's mode: whenever the weights are uniform, which they are right after every resampling, at construction and after a `weights_reset`, it returns `particles()[0]` (see [Limitations](#limitations)).

### particles

```cpp
const std::array<state_vector_t, NP>& particles() const;
```

Returns a const reference to the particle array for inspection or visualization.

### health

```cpp
pf_health health() const;
```

Returns the persistent health status. The enumerators are ordered by severity and the status latches, never returning to a lower rung on its own.

| Enumerator | Meaning |
|---|---|
| `pf_health::ok` | Every update so far conditioned the particles on its measurement |
| `pf_health::weights_reset` | A finite measurement that no particle explains reset the weights to uniform; the step succeeded with an uninformed posterior |
| `pf_health::resampling_rolled_back` | A resampling whose roughened particles would have been non-finite was undone; the step succeeded with its weighted posterior over the set before the resampling |
| `pf_health::non_finite_estimate` | A step found the particles already non-finite, which is how a poisoned control or dynamics model becomes visible |

An update refused with `non_finite_state` sets `non_finite_estimate`; the other two refusals do not set it.

## Supporting Types

### resampling_strategy

```cpp
template <typename R, typename Rng, std::size_t NP>
concept resampling_strategy = requires(const R& r, const std::array<double, NP>& weights,
                                       std::array<std::size_t, NP>& indices, Rng& rng) {
    r.resample(weights, indices, rng);
};
```

A resampler fills `indices` with NP draws from the normalized `weights`. The concept is checked with `double` weights, while the filter calls `resample` with `std::array<Scalar, NP>`, so a resampler for a `float` filter must accept both; the two built-in strategies are templates on the weight type and do.

### systematic_resampling

Low-variance resampling: one offset `u ~ U[0, 1/NP)`, then the thresholds `u + i/NP` are located in the cumulative weights by a single monotone scan. O(NP) complexity.

Header: `#include <ctrlpp/estimation/resampling/systematic_resampling.h>`

### multinomial_resampling

Multinomial resampling: NP independent uniform draws on `[0, 1)`, sorted and located in the cumulative weights. Higher variance than systematic and O(NP log NP) for the sort. Its header is not included by `particle_filter.h`; pass an instance to the second `create` overload.

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

The measurement-likelihood policy is injectable through the last class template parameter and the trailing factory argument. A policy is a callable that returns a scalar log-likelihood for a single particle's innovation:

```cpp
Scalar operator()(const Vector<Scalar, NY>& z, const Vector<Scalar, NY>& z_pred,
                  const Matrix<Scalar, NY, NY>& R_inv, Scalar log_det_2piR) const;
```

`R_inv` is the precomputed measurement precision and `log_det_2piR` is the precomputed log normalizer `log((2*pi)^NY det R)`, both derived once from the configured `R`. The default policy `detail::gaussian_likelihood<Scalar, NY>` reproduces the Gaussian log-likelihood exactly:

```cpp
Vector<Scalar, NY> innov = z - z_pred;
Scalar mahal = (innov.transpose() * R_inv * innov)(0, 0);
return Scalar{-0.5} * mahal - Scalar{0.5} * log_det_2piR;
```

A custom policy is the correct place to handle measurement topology. For a bearing sensor, for example, wrap the scalar innovation across the +/-pi cut with `std::remainder(z(0) - z_pred(0), 2*pi)` before forming the Gaussian, so a measurement near +pi and a prediction near -pi are treated as close rather than maximally distant; the usage example below does this. The policy must satisfy `pf_likelihood_model<Likelihood, Scalar, NY>`: callable as `l(z, z_pred, R_inv, log_det_2piR)` on a `const` instance, returning something convertible to `Scalar`. The value is used as a log-likelihood; `-inf` is a legitimate zero likelihood, while NaN or `+inf` refuses the update with `non_finite_result`.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'true y', '' using 1:3 with lines title 'estimated y'"

#include <ctrlpp/estimation/particle_filter.h>

#include <Eigen/Dense>

#include <cmath>
#include <random>
#include <numbers>
#include <iostream>

struct wrapped_bearing_likelihood
{
    double operator()(const ctrlpp::Vector<double, 1>& z, const ctrlpp::Vector<double, 1>& z_pred,
                      const ctrlpp::Matrix<double, 1, 1>& R_inv, double log_det_2piR) const
    {
        double innov = std::remainder(z(0) - z_pred(0), 2.0 * std::numbers::pi);
        return -0.5 * innov * R_inv(0, 0) * innov - 0.5 * log_det_2piR;
    }
};

int main()
{
    using Scalar = double;
    constexpr std::size_t NX = 2, NU = 2, NY = 1, NP = 500;
    constexpr Scalar dt = 0.1;

    // Planar position driven by a known velocity
    auto dynamics = [dt](const ctrlpp::Vector<Scalar, NX>& x,
                         const ctrlpp::Vector<Scalar, NU>& u) -> ctrlpp::Vector<Scalar, NX> {
        return x + dt * u;
    };

    // Bearing from a sensor at the origin, in (-pi, pi]
    auto measurement = [](const ctrlpp::Vector<Scalar, NX>& x) -> ctrlpp::Vector<Scalar, NY> {
        ctrlpp::Vector<Scalar, NY> z;
        z(0) = std::atan2(x(1), x(0));
        return z;
    };

    ctrlpp::pf_config<Scalar, NX, NU, NY> cfg{
        .Q = Eigen::Matrix2d::Identity() * 0.001,
        .R = (Eigen::Matrix<Scalar, 1, 1>() << 0.01).finished(),
        .x0 = Eigen::Vector2d{-3.0, 0.5},
        .P0 = Eigen::Matrix2d::Identity() * 0.5
    };

    auto created = ctrlpp::make_particle_filter<NP>(dynamics, measurement, cfg, std::mt19937_64{42},
                                                    wrapped_bearing_likelihood{});
    if (!created.has_value()) {
        std::cerr << "invalid particle filter configuration\n";
        return 1;
    }
    auto& pf = *created;

    std::mt19937 rng(123);
    std::normal_distribution<> noise(0.0, 0.1);

    ctrlpp::Vector<Scalar, NX> x_true{-3.0, 1.0};
    ctrlpp::Vector<Scalar, NU> u{0.0, -0.2};

    for (int k = 0; k < 100; ++k) {
        x_true = dynamics(x_true, u);

        ctrlpp::Vector<Scalar, NY> z;
        z(0) = std::remainder(std::atan2(x_true(1), x_true(0)) + noise(rng), 2.0 * std::numbers::pi);

        pf.predict(u);
        if (const auto stepped = pf.update(z); !stepped.has_value()) {
            std::cerr << "measurement refused at step " << k << "\n";
            continue;
        }

        auto est = pf.state();
        std::cout << k * dt << "," << x_true(1) << "," << est(1) << "\n";
    }
}
```

## Limitations

### The default likelihood assumes additive Gaussian measurement noise

The default policy `detail::gaussian_likelihood` scores each particle with the Gaussian log-density of `z - h(x)` with covariance `R`, the innovation formed as a plain vector difference. A measurement whose noise is not additive Gaussian, or that lives on a wrapped domain such as a bearing on `(-pi, pi]`, is scored wrongly: across the cut a prediction near `-pi` and a measurement near `+pi` differ by almost `2 pi` and receive almost no weight. Such measurements need a policy satisfying `pf_likelihood_model`, injected through the `Likelihood` template parameter and the trailing factory argument, as the usage example and `examples/pf/ctrlpp_pf_01_bearing_only.cpp` do. The process noise is always additive Gaussian, `x = f(x, u) + w` with `w ~ N(0, Q)`; a non-Gaussian or non-additive process noise is not supported. Measured on the usage example's track, which crosses the cut (built against libstdc++, whose normal distribution fixes the draws): with the wrapped policy the estimated `y` follows the true one with an RMS error of about 0.1; with the default likelihood it ends near `+2.9` against a true `-1.0`.

### `map_estimate()` returns the first particle after a resample

`map_estimate()` returns the particle with the largest weight, the lowest index on a tie. A resampling resets every weight to the same `1 / NP` (or `-log NP` in log mode), so until the next `update` reweights the set, `map_estimate()` is `particles()[0]`: one roughened copy chosen by its position in the array, not the mode of the posterior. The same holds at construction and after a `weights_reset`. With `extraction` set to `map`, `state()` returns that particle on every step that resampled. No kernel density or other mode estimate is computed; the weighted mean, or a mode formed by the caller from `particles()`, is the alternative.

## See Also

- [ukf](ukf.md)<br/> sigma-point filter for moderate nonlinearity
- [ekf](ekf.md)<br/> linearization-based filter for smooth systems
- [observer-policy](observer-policy.md)<br/> concept satisfied by this type
- [background/particle-filter](../../background/particle-filter.md)<br/> SIR derivation
- [guides/estimation/observer-controller](../../guides/estimation/observer-controller.md)<br/> composing observers with controllers
