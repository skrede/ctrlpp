# Your First Estimator

This tutorial builds a Kalman filter that tracks a constant-velocity target
using noisy position-only measurements. The predict-update cycle is the
foundation of all estimation in ctrlpp.

**Prerequisites:** ctrlpp installed per [Getting Started](../../getting-started.md).

## The System

A constant-velocity model in one dimension:

```
state  = [position, velocity]
input  = [] (no control input)
output = [position] (we only measure position)
```

The discrete state-space matrices are:

```
A = [1  dt]    B = [0]    C = [1  0]
    [0  1 ]        [0]
```

## Complete Program

```cpp
// Usage: ./your_first_estimator | gnuplot -p -e "set datafile separator ','; plot '-' skip 1 using 1:2 with lines title 'true', '' using 1:4 with lines title 'estimated', '' using 1:6 with points pt 7 ps 0.3 title 'measured'"
#include <ctrlpp/estimation/kalman.h>
#include <ctrlpp/model/state_space.h>

#include <Eigen/Dense>

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <random>

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;
    constexpr std::size_t NY = 1;
    constexpr double dt = 0.1;
    constexpr std::size_t n_steps = 200;

    // Discrete state-space: constant velocity model
    ctrlpp::discrete_state_space<double, NX, NU, NY> sys{
        .A = (Eigen::Matrix2d() << 1.0, dt, 0.0, 1.0).finished(),
        .B = Eigen::Vector2d::Zero(),
        .C = (Eigen::Matrix<double, 1, 2>() << 1.0, 0.0).finished(),
        .D = Eigen::Matrix<double, 1, 1>::Zero()};

    // Process and measurement noise covariances
    Eigen::Matrix2d Q = Eigen::Matrix2d::Identity() * 0.01;
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    // Initial state estimate and covariance
    Eigen::Vector2d x0 = Eigen::Vector2d::Zero();
    Eigen::Matrix2d P0 = Eigen::Matrix2d::Identity() * 10.0;

    ctrlpp::kalman_filter<double, NX, NU, NY> kf(
        sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0});

    // True state: position=0, velocity=1 (constant velocity)
    Eigen::Vector2d x_true;
    x_true << 0.0, 1.0;

    std::mt19937 gen(42);
    std::normal_distribution<double> meas_noise(0.0, 1.0);

    std::cout << "time,true_pos,true_vel,est_pos,est_vel,meas_pos\n";

    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();

    for (std::size_t k = 0; k < n_steps; ++k)
    {
        double t = static_cast<double>(k) * dt;

        // Predict
        kf.predict(u);

        // Propagate true state
        x_true = sys.A * x_true;

        // Noisy measurement
        double z_val = x_true(0) + meas_noise(gen);
        Eigen::Matrix<double, 1, 1> z;
        z << z_val;

        // Update
        kf.update(z);

        auto est = kf.state();

        std::cout << std::fixed << std::setprecision(4)
                  << t << ","
                  << x_true(0) << "," << x_true(1) << ","
                  << est(0) << "," << est(1) << ","
                  << z_val << "\n";
    }
}
```

## The Predict-Update Cycle

Every step in the loop follows the same two-phase pattern:

1. **Predict** -- `kf.predict(u)` propagates the state estimate forward using
   the system model. The covariance grows (uncertainty increases).

2. **Update** -- `kf.update(z)` incorporates the new measurement. The Kalman
   gain balances the predicted estimate against the measurement based on their
   relative uncertainties. The covariance shrinks.

This cycle repeats at every time step. The filter converges quickly -- within
a few steps the estimated velocity (which is not directly measured) tracks the
true velocity closely.

## Tuning the Noise Covariances

- **Q** (process noise) -- increase if the model is less trustworthy.
  Larger Q makes the filter respond faster to changes but be noisier.
- **R** (measurement noise) -- increase if measurements are noisy. Larger R
  makes the filter smoother but slower to react.

## Next Steps

- [Kalman Filter API](../../api/estimation/kalman.md) -- full interface reference
- [EKF API](../../api/estimation/ekf.md) -- for nonlinear systems
- [Observer-Controller Composition](../estimation/observer-controller.md) --
  feeding estimates into a controller
- [Kalman Theory](../../background/kalman.md) -- derivation and
  optimality properties
