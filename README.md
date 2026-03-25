# ctrlpp
[![Linux](https://github.com/skrede/ctrlpp/actions/workflows/linux.yml/badge.svg?branch=master)](https://github.com/skrede/ctrlpp/actions/workflows/linux.yml)
[![macOS](https://github.com/skrede/ctrlpp/actions/workflows/macos.yml/badge.svg?branch=master)](https://github.com/skrede/ctrlpp/actions/workflows/macos.yml)
[![Windows](https://github.com/skrede/ctrlpp/actions/workflows/windows.yml/badge.svg?branch=master)](https://github.com/skrede/ctrlpp/actions/workflows/windows.yml)
[![codecov](https://codecov.io/gh/skrede/ctrlpp/branch/master/graph/badge.svg)](https://codecov.io/gh/skrede/ctrlpp)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue)](LICENSE)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://en.cppreference.com/w/cpp/20)

**ctrlpp** is a C++20 control systems library with policy-based composition and concept-constrained interfaces. Header-only, Eigen-backed, and designed for hot-path use in real-time systems. PID controllers compose from orthogonal policies (anti-windup, derivative filtering, rate limiting); estimators and MPC/MHE inject solver backends through concepts; system identification runs online or offline with unified result types.

**NB:** This library is still under development and have not yet been tested extensively in real-world scenarios -- beyond the (substantial) test suite under `tests/`. Reports and experiences use and testing of this library will be appreciated.

## Features

- **Policy-based PID** -- compose anti-windup, derivative filtering, setpoint filtering, velocity form, ISA form, feed-forward, and rate limiting from orthogonal policy types.
- **Estimation** -- Kalman, Luenberger, EKF, UKF, particle filter, MEKF, manifold UKF, and complementary filter with a unified observer concept interface.
- **Model predictive control** -- linear MPC (OSQP) and nonlinear MPC (NLopt) with terminal constraints, soft constraints, and delta-u limiting.
- **Moving horizon estimation** -- linear MHE (OSQP) and nonlinear MHE (NLopt) with arrival cost and box constraints.
- **Signal processing** -- biquad IIR sections (Butterworth, Chebyshev), FIR filters, and cascaded filter chains.
- **System identification** -- RLS, batch/recursive ARX, and N4SID subspace identification with fit metrics.
- **Lie group utilities** -- SO(3) quaternion exponential/logarithm maps for attitude estimation.
- **Model utilities** -- state-space and transfer function representations, discretisation, conversion, stability analysis, and C++20 concepts for dynamics, measurement, and constraint models.

## Quick Start

```cpp
#include "ctrlpp/control/pid.h"

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
    constexpr double setpoint = 1.0;
    constexpr double dt = 0.01;

    for (double t = 0.0; t < 10.0; t += dt)
    {
        auto u = ctrl.compute(Vec::Constant(setpoint), Vec::Constant(y), dt);
        y = 0.9 * y + 0.1 * u[0];
    }

    std::cout << "final output: " << y << "\n";
}
```

`pid<double, 1, 1, 1>` is a SISO PID with double precision. `config_type` sets gains and output limits. `compute()` takes setpoint, measurement, and timestep, returning the control signal. Add policies for richer behaviour:

```cpp
#include "ctrlpp/control/pid.h"

using Pid = ctrlpp::pid<double, 1, 1, 1,
    ctrlpp::anti_windup_policy,
    ctrlpp::derivative_filter_policy,
    ctrlpp::rate_limit_policy>;
```

## CMake Integration

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

### find_package

```cmake
find_package(ctrlpp CONFIG REQUIRED)
target_link_libraries(my_app PRIVATE ctrlpp::ctrlpp)
```

### Optional solver backends

OSQP and NLopt are fetched automatically when MPC/MHE headers are used. To disable:

```cmake
set(CTRLPP_ENABLE_OSQP OFF)
set(CTRLPP_ENABLE_NLOPT OFF)
```

## Documentation

- [Getting Started](docs/getting-started.md) -- Install ctrlpp and run your first PID controller
- [Guides](docs/guides/README.md) -- Tutorials and deep dives
- [API Reference](docs/README.md#api-reference) -- Full type documentation
- [Reference](docs/reference/README.md) -- Theory and mathematical background

## License

Apache 2.0 License -- see [LICENSE](LICENSE) for details.
