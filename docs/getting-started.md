# Getting Started

## Prerequisites

- C++23 compiler: GCC 14+, Clang 18+, MSVC 17.10+
- CMake 3.25+
- Eigen 3.4+ (fetched automatically via FetchContent)

## Installation

### FetchContent (recommended)

```cmake
include(FetchContent)
FetchContent_Declare(
    ctrlpp
    GIT_REPOSITORY https://github.com/skrede/ctrlpp.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(ctrlpp)

target_link_libraries(my_app PRIVATE ctrlpp::ctrlpp)
```

This pulls ctrlpp and its Eigen dependency automatically. No manual
installation required.

### find_package

```cmake
find_package(ctrlpp CONFIG REQUIRED)
target_link_libraries(my_app PRIVATE ctrlpp::ctrlpp)
```

## Your first PID controller

The following program creates a SISO PID controller and runs it against a
simple first-order plant for 10 seconds, printing the response every half
second.

```cpp
#include "ctrlpp/pid.h"

#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1, 1, 1>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.0);
    cfg.ki = Vec::Constant(1.0);
    cfg.kd = Vec::Constant(0.0);
    cfg.output_min = Vec::Constant(-10.0);
    cfg.output_max = Vec::Constant(10.0);

    Pid ctrl(cfg);

    double y = 0.0;
    constexpr double dt = 0.01;

    for (int i = 0; i < 1000; ++i)
    {
        auto u = ctrl.compute(Vec::Constant(1.0), Vec::Constant(y), dt);
        y = 0.9 * y + 0.1 * u(0);

        if (i % 50 == 0)
            std::cout << i * dt << "s: y=" << y << " u=" << u(0) << "\n";
    }
}
```

`pid<double, 1, 1, 1>` is a single-input single-output PID controller with
`double` precision. The three `1`s are the state, input, and output dimensions
-- all one for SISO. For MIMO systems, increase these to match your plant.

`config_type` holds the controller gains and output limits. `kp`, `ki`, and
`kd` are Eigen vectors sized to the output dimension. `output_min` and
`output_max` clamp the control signal.

`compute()` takes the setpoint, the current measurement, and the timestep as
arguments. It returns the control signal as an Eigen vector. The controller
tracks the integral and derivative state internally between calls.

## What's next

- [Your First Estimator](guides/intro/your-first-estimator.md) -- add an observer to your control loop
- [Your First MPC](guides/intro/your-first-mpc.md) -- model predictive control with constraints
- [PID Composition](guides/pid/composition.md) -- add anti-windup, derivative filtering, and more
- [API Reference](README.md#api-reference) -- full type documentation
