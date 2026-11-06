# observer_policy

Concept definitions for observer types in ctrlpp. The `ObserverPolicy` concept states the minimal interface every ctrlpp observer has, so user code can write one loop for any of them; no library component is constrained on it, and controllers take a state, not an observer. The `CovarianceObserver` refinement adds covariance and innovation access. A `null_observer` type satisfies `ObserverPolicy` with no-op operations, for generic user code that is instantiated without an observer.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::ObserverPolicy` | `#include <ctrlpp/estimation/observer_policy.h>` |
| `ctrlpp::CovarianceObserver` | `#include <ctrlpp/estimation/observer_policy.h>` |
| `ctrlpp::null_observer` | `#include <ctrlpp/estimation/observer_policy.h>` |
| All of the above | `#include <ctrlpp/observer_policy.h>` (convenience) |

## Concepts

### ObserverPolicy

```cpp
template <typename O>
concept ObserverPolicy = requires {
    typename O::observer_tag;
    typename O::state_vector_t;
    typename O::input_vector_t;
    typename O::output_vector_t;
} && requires(O obs, const typename O::input_vector_t& u, const typename O::output_vector_t& z) {
    obs.predict(u);
    obs.update(z);
    { obs.state() } -> std::convertible_to<const typename O::state_vector_t&>;
};
```

Any type satisfying `ObserverPolicy` provides:

| Requirement | Description |
|-------------|-------------|
| `observer_tag` | A member type; the library's observers alias a distinct struct per type and `null_observer` aliases `void`. Nothing in the library reads it |
| `state_vector_t` | The state vector type |
| `input_vector_t` | The control input type |
| `output_vector_t` | The measurement type |
| `predict(u)` | Propagate state estimate forward one step; any return type |
| `update(z)` | Incorporate a measurement; any return type (the library's observers return `ctrlpp::expected<void, E>`) |
| `state()` | Return the current state estimate |

### CovarianceObserver

```cpp
template <typename O>
concept CovarianceObserver = ObserverPolicy<O> && requires(const O& obs) {
    { obs.covariance() };
    { obs.innovation() };
};
```

Refines `ObserverPolicy` with covariance and innovation access; the concept requires the two calls to compile and constrains neither return type. Satisfied by `kalman_filter`, `ekf`, `ukf`, `mekf`, `manifold_ukf`, `mhe` and `nmhe`.

## Types Satisfying ObserverPolicy

| Type | CovarianceObserver | Header |
|------|-------------------|--------|
| `kalman_filter` | yes | `estimation/kalman.h` |
| `luenberger_observer` | no | `estimation/luenberger.h` |
| `ekf` | yes | `estimation/ekf.h` |
| `ukf` | yes | `estimation/ukf.h` |
| `particle_filter` | no (see below) | `estimation/particle_filter.h` |
| `mekf` | yes | `estimation/mekf.h` |
| `manifold_ukf` | yes | `estimation/manifold_ukf.h` |
| `complementary_filter` | no | `estimation/complementary_filter.h` |
| `mhe` | yes | `mhe.h` |
| `nmhe` | yes | `nmhe.h` |
| `null_observer` | no | `estimation/observer_policy.h` |

`particle_filter` does expose a `covariance()` accessor returning the weighted posterior covariance of its particle set, but it carries no innovation, so it does not satisfy the refinement. The concept requires both.

## null_observer

A type with no-op operations and `std::monostate` as every vector type; `state()` returns a reference to a function-local static `std::monostate`. It satisfies `ObserverPolicy` and not `CovarianceObserver`.

```cpp
struct null_observer
{
    using observer_tag    = void;
    using state_vector_t  = std::monostate;
    using input_vector_t  = std::monostate;
    using output_vector_t = std::monostate;

    void predict(const std::monostate&) {}
    void update(const std::monostate&) {}
    const std::monostate& state() const;
};
```

## Usage Example

```cpp
#include <ctrlpp/estimation/observer_policy.h>
#include <ctrlpp/estimation/kalman.h>
#include <ctrlpp/control/lqr.h>
#include <ctrlpp/model/state_space.h>

#include <Eigen/Dense>

#include <iostream>

// Generic function accepting any observer
template <ctrlpp::ObserverPolicy Observer>
void run_loop(Observer& obs, auto& ctrl, auto& sys, int steps)
{
    using input_t = typename Observer::input_vector_t;
    using output_t = typename Observer::output_vector_t;

    Eigen::Vector2d x_true;
    x_true << 1.0, 0.0;

    for (int k = 0; k < steps; ++k) {
        input_t u = ctrl.compute(obs.state());

        x_true = sys.A * x_true + sys.B * u;
        output_t z = sys.C * x_true;

        obs.predict(u);
        // The concept does not constrain update's return type, so a generic
        // observer loop that wants the diagnosis constrains it itself; the
        // library's own observers all return ctrlpp::expected<void, E>.
        if(!obs.update(z))
        {
            std::cerr << "observer rejected the measurement at step " << k << "\n";
            return;
        }

        std::cout << k << "," << x_true(0) << "," << obs.state()(0) << "\n";
    }
}

int main()
{
    using Scalar = double;
    constexpr std::size_t NX = 2, NU = 1, NY = 1;

    ctrlpp::discrete_state_space<Scalar, NX, NU, NY> sys{};
    sys.A << 1.0, 0.1, 0.0, 1.0;
    sys.B << 0.005, 0.1;
    sys.C << 1.0, 0.0;
    sys.D << 0.0;

    auto K_opt = ctrlpp::lqr_gain<Scalar, NX, NU>(
        sys.A, sys.B,
        Eigen::Matrix2d::Identity() * 10.0,
        Eigen::Matrix<Scalar, 1, 1>::Identity());
    if (!K_opt.has_value())
        return 1;

    ctrlpp::lqr<Scalar, NX, NU> ctrl(*K_opt);

    auto kf_result = ctrlpp::kalman_filter<Scalar, NX, NU, NY>::create(sys, {});
    if (!kf_result)
        return 1;
    auto& kf = *kf_result;

    run_loop(kf, ctrl, sys, 100);
    return 0;
}
```

## See Also

- [kalman](kalman.md)<br/> linear Kalman filter
- [ekf](ekf.md)<br/> extended Kalman filter
- [ukf](ukf.md)<br/> unscented Kalman filter
- [particle-filter](particle-filter.md)<br/> particle filter
- [luenberger](luenberger.md)<br/> Luenberger observer
- [guides/estimation/observer-controller](../../guides/estimation/observer-controller.md)<br/> composition patterns
