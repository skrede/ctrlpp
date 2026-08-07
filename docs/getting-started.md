# Getting Started

## Prerequisites

- C++20 compiler: GCC 13+, Clang 18+, MSVC 17.10+, Xcode 15.4+<br/>
  Fallible solvers return `ctrlpp::expected`, which resolves to `std::expected`
  automatically when compiled at C++23 or later and to an in-library C++20
  fallback with the same call surface otherwise.
- CMake 3.28+ to build ctrlpp<br/>
  That floor is what building ctrlpp itself requires. A project that only consumes an
  installed ctrlpp may declare a lower `cmake_minimum_required`; the installed package
  imposes no such floor on it.
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

## Versioning and compatibility

An install answers a versioned request. An install of `0.3.6` satisfies a request for the same
major and minor version with an equal or lower patch, and refuses a request for any other minor
version. So a request for `0.3.0` succeeds against it:

```cmake
find_package(ctrlpp 0.3.0 CONFIG REQUIRED)
target_link_libraries(my_app PRIVATE ctrlpp::ctrlpp)
```

while `find_package(ctrlpp 0.1.0 CONFIG REQUIRED)` against the same install fails with `Could not
find a configuration file for package "ctrlpp" that is compatible with requested version "0.1.0"`.
Before 1.0 that is the honest contract, because a minor bump is free to break the API. Past 1.0 it
is stricter than semantic versioning requires, since a 1.1 install would refuse a request for 1.0
even though nothing broke; the intended replacement then is CMake's semantic-version compatibility
mode, which is not used today because generating it requires CMake 4.4 and this project asks for
3.28.

## Solver backend components

Backends are requested as package components. The component names are `osqp`, `nlopt` and `argmin`,
and each one brings a target of the same name under the `ctrlpp::` namespace:

```cmake
find_package(ctrlpp CONFIG REQUIRED COMPONENTS nlopt)
target_link_libraries(my_app PRIVATE ctrlpp::nlopt)
```

`ctrlpp::osqp` and `ctrlpp::argmin` are requested the same way, naming `osqp` or `argmin` as the
component. An install answers a component request only if it was built with that backend; asking
for one it does not carry fails at `find_package` time and names the component:

```
this ctrlpp install was built without the requested component(s): osqp.
Rebuild ctrlpp with the matching CTRLPP_BUILD_<COMPONENT> option on against
a discoverable dependency and reinstall.
```

The installed package resolves each backend's own dependency on your behalf, so your environment
has to make that dependency discoverable: an OSQP installation for `osqp`, NLopt 2.10 or newer for
`nlopt`, and Argmin 0.3 or newer for `argmin`. Two of those requirements are looser than they read.
The OSQP one carries no version floor, because OSQP's package publishes no project version and every
install of it reports `0.0.0`. The Argmin floor is satisfied by any 0.x install today, because
Argmin's package declares major-version compatibility; it starts to bite once Argmin narrows that
mode.

## Your first PID controller

The following program creates a SISO PID controller and runs it against a
simple first-order plant for 10 seconds, printing the response every half
second.

```cpp
// Usage: ./first_pid | gnuplot -p -e "set datafile separator ' '; plot '-' using 1:2 with lines title 'output'"
#include <ctrlpp/control/pid.h>

#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1>;
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
        // A refused cycle produced no command; the caller decides what the
        // actuator does. This sample stops.
        auto step = ctrl.compute(Vec::Constant(1.0), Vec::Constant(y), dt);
        if (!step.has_value()) return 1;
        const auto& u = *step;
        y = 0.9 * y + 0.1 * u(0);

        if (i % 50 == 0)
            std::cout << i * dt << "s: y=" << y << " u=" << u(0) << "\n";
    }
}
```

`pid<double, 1>` is a single-input single-output PID controller with
`double` precision. The three `1`s are the state, input, and output dimensions
-- all one for SISO. For MIMO systems, increase these to match your plant.

`config_type` holds the controller gains and output limits. `kp`, `ki`, and
`kd` are Eigen vectors sized to the output dimension. `output_min` and
`output_max` clamp the control signal.

`compute()` takes the setpoint, the current measurement, and the timestep as
arguments. It returns the control signal as an Eigen vector. The controller
tracks the integral and derivative state internally between calls.

## What's next

- [Your First Estimator](guides/intro/your-first-estimator.md)<br/> add an observer to your control loop
- [Your First MPC](guides/intro/your-first-mpc.md)<br/> model predictive control with constraints
- [PID Composition](guides/pid/composition.md)<br/> add anti-windup, derivative filtering, and more
- [API Reference](README.md#api-reference)<br/> full type documentation
