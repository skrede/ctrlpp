# observer_policy

Concept definitions for observer types in ctrlpp. The `ObserverPolicy` concept defines the minimal interface that any state observer must satisfy for composition with controllers and other library components. The `CovarianceObserver` refinement adds covariance and innovation access for types that maintain probabilistic state estimates. A `null_observer` sentinel type is also provided for controllers that do not need an observer.

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
| `observer_tag` | A unique tag type for type identification |
| `state_vector_t` | The state vector type |
| `input_vector_t` | The control input type |
| `output_vector_t` | The measurement type |
| `predict(u)` | Propagate state estimate forward one step |
| `update(z)` | Incorporate a measurement |
| `state()` | Return the current state estimate |

### CovarianceObserver

```cpp
template <typename O>
concept CovarianceObserver = ObserverPolicy<O> && requires(const O& obs) {
    { obs.covariance() };
    { obs.innovation() };
};
```

Refines `ObserverPolicy` with covariance and innovation access. Satisfied by `kalman_filter`, `ekf`, `ukf`, `mekf`, and `manifold_ukf`.

## Types Satisfying ObserverPolicy

| Type | CovarianceObserver | Header |
|------|-------------------|--------|
| `kalman_filter` | yes | `estimation/kalman.h` |
| `luenberger_observer` | no | `estimation/luenberger.h` |
| `ekf` | yes | `estimation/ekf.h` |
| `ukf` | yes | `estimation/ukf.h` |
| `particle_filter` | no | `estimation/particle_filter.h` |
| `mekf` | yes | `estimation/mekf.h` |
| `manifold_ukf` | yes | `estimation/manifold_ukf.h` |
| `complementary_filter` | no | `estimation/complementary_filter.h` |
| `null_observer` | no | `estimation/observer_policy.h` |

## null_observer

A sentinel type for controllers that do not use an observer. All operations are no-ops with `std::monostate` as the vector type.

```cpp
struct null_observer
{
    using observer_tag    = void;
    using state_vector_t  = std::monostate;
    using input_vector_t  = std::monostate;
    using output_vector_t = std::monostate;

    void predict(const std::monostate&) {}
    void update(const std::monostate&) {}
    auto state() const -> const std::monostate&;
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
        auto x_est = obs.state();
        auto u = ctrl.compute(x_est);

        output_t z = sys.C * x_true;
        x_true = sys.A * x_true + sys.B * u;

        obs.predict(u);
        obs.update(z);

        std::cout << k << "," << x_true(0) << "," << x_est(0) << "\n";
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

    ctrlpp::lqr<Scalar, NX, NU> ctrl(*K_opt);

    ctrlpp::kalman_filter<Scalar, NX, NU, NY> kf(sys, {});

    run_loop(kf, ctrl, sys, 100);
}
```

## See Also

- [kalman](kalman.md) -- linear Kalman filter
- [ekf](ekf.md) -- extended Kalman filter
- [ukf](ukf.md) -- unscented Kalman filter
- [particle-filter](particle-filter.md) -- particle filter
- [luenberger](luenberger.md) -- Luenberger observer
- [guides/estimation/observer-controller](../../guides/estimation/observer-controller.md) -- composition patterns
