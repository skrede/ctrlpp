# Observer-Controller Composition

In many real systems, not all states are directly measurable. An observer
estimates the full state vector from partial measurements, and a controller
uses that estimate. ctrlpp keeps observers and controllers as separate objects
that compose through the `ObserverPolicy` concept -- no coupling, no
inheritance.

## The Pattern

```
measurement --> [Observer] --> estimated_state --> [Controller] --> control
                    ^                                                  |
                    |-------------- control (for predict) -------------|
```

The observer's `state()` feeds the controller's `compute()` or `solve()`.
The controller's output feeds back into the observer's `predict()` on the
next step. This is standard separation of estimation and control.

## Kalman Filter + LQR

This example uses a Kalman filter to estimate the full state of a cart-pendulum
from position and angle measurements, then feeds the estimate into an LQR
controller.

```cpp
#include "ctrlpp/estimation/kalman.h"
#include "ctrlpp/model/discretise.h"
#include "ctrlpp/control/lqr.h"
#include "ctrlpp/model/propagate.h"
#include "ctrlpp/model/state_space.h"

#include <iomanip>
#include <iostream>

int main()
{
    using Scalar = double;
    constexpr std::size_t NX = 4;
    constexpr std::size_t NU = 1;
    constexpr std::size_t NY = 2;

    constexpr Scalar M_cart = 1.0;
    constexpr Scalar m_pend = 0.1;
    constexpr Scalar l = 0.5;
    constexpr Scalar g = 9.81;
    constexpr Scalar dt = 0.05;
    constexpr Scalar duration = 20.0;

    constexpr Scalar denom = M_cart + m_pend;

    // Linearised cart-pendulum continuous model
    ctrlpp::continuous_state_space<Scalar, NX, NU, NY> sys_c{};
    sys_c.A << 0.0, 1.0, 0.0, 0.0,
               0.0, 0.0, -m_pend * g / denom, 0.0,
               0.0, 0.0, 0.0, 1.0,
               0.0, 0.0, (denom * g) / (denom * l), 0.0;
    sys_c.B << 0.0, 1.0 / denom, 0.0, -1.0 / (denom * l);
    sys_c.C << 1.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 1.0, 0.0;
    sys_c.D.setZero();

    auto sys_d = ctrlpp::discretise(ctrlpp::zoh{}, sys_c, dt);

    // LQR gain
    Eigen::Matrix<Scalar, 4, 4> Q_lqr = Eigen::Matrix<Scalar, 4, 4>::Zero();
    Q_lqr(0, 0) = 10.0;
    Q_lqr(1, 1) = 1.0;
    Q_lqr(2, 2) = 100.0;
    Q_lqr(3, 3) = 10.0;
    Eigen::Matrix<Scalar, 1, 1> R_lqr;
    R_lqr << 1.0;

    auto K_opt = ctrlpp::lqr_gain<Scalar, NX, NU>(sys_d.A, sys_d.B, Q_lqr, R_lqr);
    ctrlpp::lqr<Scalar, NX, NU> controller(*K_opt);

    // Kalman filter
    Eigen::Matrix<Scalar, 4, 4> Q_proc = Eigen::Matrix<Scalar, 4, 4>::Identity() * 0.01;
    Eigen::Matrix<Scalar, 2, 2> R_meas = Eigen::Matrix<Scalar, 2, 2>::Identity() * 0.01;
    Eigen::Matrix<Scalar, 4, 1> x0_est = Eigen::Matrix<Scalar, 4, 1>::Zero();
    Eigen::Matrix<Scalar, 4, 4> P0 = Eigen::Matrix<Scalar, 4, 4>::Identity();

    ctrlpp::kalman_filter<Scalar, NX, NU, NY> kf(
        sys_d, {.Q = Q_proc, .R = R_meas, .x0 = x0_est, .P0 = P0});

    // True initial state (unknown to the observer)
    Eigen::Matrix<Scalar, 4, 1> x_true;
    x_true << 0.1, 0.0, 0.05, 0.0;

    std::cout << "time,true_pos,est_pos,true_angle,est_angle,control\n";

    for (Scalar t = 0.0; t < duration; t += dt)
    {
        // Observer estimate feeds into controller
        auto x_est = kf.state();
        auto u = controller.compute(x_est);

        // Measurement (position + angle)
        Eigen::Matrix<Scalar, 2, 1> z = sys_d.C * x_true;

        std::cout << std::fixed << std::setprecision(4)
                  << t << "," << x_true(0) << "," << x_est(0) << ","
                  << x_true(2) << "," << x_est(2) << ","
                  << u(0) << "\n";

        // Predict observer with applied control
        kf.predict(u);

        // Propagate true state
        x_true = ctrlpp::propagate(sys_d, x_true, u);

        // Update observer with measurement
        kf.update(z);
    }
}
```

## Kalman Filter + MPC

The same pattern works with MPC. The Kalman filter provides the estimated state
to `controller.solve()`:

```cpp
// In the control loop:
kf.predict(u);
x_true = ctrlpp::propagate(mpc_sys, x_true, u);
Eigen::Matrix<double, 1, 1> z = obs_sys.C * x_true;
kf.update(z);

auto u_opt = controller.solve(kf.state());
```

See the [MPC tutorial](../intro/your-first-mpc.md) for the full MPC setup.

## The ObserverPolicy Concept

Any type satisfying `ObserverPolicy` can serve as an observer:

```cpp
template <typename O>
concept ObserverPolicy = requires {
    typename O::observer_tag;
    typename O::state_vector_t;
    typename O::input_vector_t;
    typename O::output_vector_t;
} && requires(O obs,
              const typename O::input_vector_t& u,
              const typename O::output_vector_t& z) {
    obs.predict(u);
    obs.update(z);
    { obs.state() } -> std::convertible_to<const typename O::state_vector_t&>;
};
```

All ctrlpp observers -- `kalman_filter`, `luenberger_observer`, `ekf`, `ukf`,
`particle_filter` -- satisfy this concept. You can swap observers without
changing the controller code.

## Which Observer to Use

| Observer             | Use case                                   |
| -------------------- | ------------------------------------------ |
| `kalman_filter`      | Linear systems, Gaussian noise             |
| `luenberger_observer`| Linear systems, no noise model needed      |
| `ekf`                | Mildly nonlinear, differentiable dynamics  |
| `ukf`                | Nonlinear, no Jacobians needed             |
| `particle_filter`    | Strongly nonlinear, non-Gaussian           |

## Next Steps

- [ObserverPolicy API](../../estimation/observer-policy.md) -- concept
  definition and null_observer
- [Kalman Filter API](../../estimation/kalman.md) -- linear observer
- [EKF API](../../estimation/ekf.md) -- nonlinear observer
- [LQR API](../../control/lqr.md) -- linear-quadratic regulator
- [Kalman Theory](../../reference/kalman-theory.md) -- optimality and
  convergence
