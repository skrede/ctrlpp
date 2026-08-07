# ctrlpp
[![Linux](https://github.com/skrede/ctrlpp/actions/workflows/linux.yml/badge.svg?branch=master)](https://github.com/skrede/ctrlpp/actions/workflows/linux.yml)
[![macOS](https://github.com/skrede/ctrlpp/actions/workflows/macos.yml/badge.svg?branch=master)](https://github.com/skrede/ctrlpp/actions/workflows/macos.yml)
[![Windows](https://github.com/skrede/ctrlpp/actions/workflows/windows.yml/badge.svg?branch=master)](https://github.com/skrede/ctrlpp/actions/workflows/windows.yml)
[![codecov](https://codecov.io/gh/skrede/ctrlpp/branch/master/graph/badge.svg)](https://codecov.io/gh/skrede/ctrlpp)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue)](LICENSE)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://en.cppreference.com/w/cpp/20)

**ctrlpp** is a C++20 control systems library with policy-based composition and concept-constrained interfaces. Header-only and Eigen-backed. PID controllers compose from orthogonal policies (anti-windup, derivative filtering, rate limiting); estimators and MPC/MHE inject solver backends through concepts; system identification runs online or offline with unified result types.

**NB:** This library is still under development and has not undergone rigorous real-world testing beyond the extensive test suite under `tests/`. Reports and experiences from use or testing of this library will be appreciated.

## Features

- **Policy-based PID**<br/> Compose anti-windup, derivative filtering, setpoint filtering, velocity form, ISA form, feed-forward, and rate limiting from orthogonal policy types.
- **Estimation**<br/> Kalman, Luenberger, EKF, UKF, particle filter, MEKF, manifold UKF, and complementary filter with a unified observer concept interface.
- **Model predictive control**<br/> linear MPC (OSQP) and nonlinear MPC (NLopt) with terminal constraints, soft constraints, and delta-u limiting.
- **Moving horizon estimation**<br/> linear MHE (OSQP) and nonlinear MHE (NLopt) with arrival cost and box constraints.
- **Signal processing**<br/> biquad IIR sections (Butterworth, Chebyshev), FIR filters, and cascaded filter chains.
- **System identification**<br/> RLS, batch/recursive ARX, and N4SID subspace identification with fit metrics.
- **Lie group utilities**<br/> SO(3) quaternion exponential/logarithm maps for attitude estimation.
- **Model utilities**<br/> state-space and transfer function representations, discretisation, conversion, stability analysis, and C++20 concepts for dynamics, measurement, and constraint models.

## Quick Start

```cpp
#include "ctrlpp/control/pid.h"

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

`pid<double, 1>` is a SISO PID with double precision. `config_type` sets gains and output limits. `compute()` takes setpoint, measurement, and timestep, returning the control signal. Add policies for richer behaviour:

```cpp
#include "ctrlpp/control/pid.h"

using Pid = ctrlpp::pid<double, 1,
    ctrlpp::anti_windup<ctrlpp::back_calc>,
    ctrlpp::deriv_filter,
    ctrlpp::rate_limit>;
```

The four leading template arguments (`Scalar, NX, NU, NY`) are the current arity; a future
release will collapse the three size parameters behind a single output-count parameter and
keep a deprecated alias for the form shown here, so existing code keeps compiling.
`anti_windup` selects an anti-windup strategy tag: back-calculation (`anti_windup<back_calc>`,
shown above), clamping (`anti_windup<clamping>`), or conditional integration
(`anti_windup<conditional_integration>`).

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

OSQP, NLopt, and Argmin are opt-in backends: each is off by default and is requested by its own
option.

```cmake
set(CTRLPP_BUILD_OSQP ON)
set(CTRLPP_BUILD_NLOPT ON)
set(CTRLPP_BUILD_ARGMIN ON)
set(CTRLPP_CMAKE_FETCH_DEPS ON)
```

Every dependency, Eigen included, is declared once and carries both acquisition routes;
`CTRLPP_CMAKE_FETCH_DEPS` selects between them. With it off (the default), a discoverable
installation of the dependency is used, and the pinned source is fetched only when that lookup finds
nothing. With it on, the pinned source is built inside this build tree and no installed copy is
consulted.

A dependency acquired by fetching cannot be shipped. A target built inside this build tree belongs
to no export set, so ctrlpp's own export cannot reference it, and asking for a fetch and an install
in the same build is a configure error naming every dependency in that state rather than a silent
downgrade:

```
ctrlpp cannot be installed.  These dependencies were fetched into this
build and install nothing, so they belong to no export set and ctrlpp's own
export cannot reference them: Eigen3, osqp.  Either supply a discoverable
version of each one and configure with CTRLPP_CMAKE_FETCH_DEPS=OFF, or
configure with CTRLPP_INSTALL=OFF if this build was never meant to ship.
```

Eigen is subject to the same rule, so a fetching build cannot be installed at all. Shipping ctrlpp,
with a backend or without one, therefore requires a discoverable installation of every dependency it
uses:

```sh
cmake -S . -B build -DCTRLPP_BUILD_OSQP=ON -DCTRLPP_CMAKE_FETCH_DEPS=OFF -DCTRLPP_INSTALL=ON
cmake --install build --prefix /your/prefix
```

`CTRLPP_INSTALL` decides whether ctrlpp reaches the install prefix at all. It defaults to on for a
top-level build and off when ctrlpp is added as a subproject, so a parent project's `cmake --install`
carries only what that parent asked for; a parent that does want ctrlpp installed sets the option
before adding it.

`CTRLPP_ARGMIN_GIT_TAG` and `CTRLPP_ARGMIN_SOURCE_DIR` pin the Argmin backend to a specific
git tag or a local source checkout, respectively.

### Exception posture

ctrlpp is **consumer-flag-agnostic**: it forces no `-fno-exceptions` / `-fno-rtti` on any
installed or interface target, and `config.h` auto-detects `__cpp_exceptions` to adapt to
whatever you compile with. Every type is built through a fallible factory returning
`ctrlpp::expected`, and no construction path is gated on exceptions. The only wrappers still
gated are the `setup(problem)` convenience overloads on the optional OSQP and NLopt backend
adapters, whose fallible `try_setup` counterparts are unconditional.

As a self-imposed compatibility guarantee, ctrlpp's **own** tests and benches dogfood the
throw-free discipline. `CTRLPP_TESTS_WITH_EXCEPTIONS` selects the build tree:

- **OFF (default, the `dev` preset)** — ctrlpp's test targets **and** Catch2 compile
  `-fno-exceptions -fno-rtti` (Catch2 with `CATCH_CONFIG_DISABLE_EXCEPTIONS`). The
  allocation-free no-malloc suite runs on this build, so zero-alloc and no-throw are proven
  together. A blocking CI job gates it.
- **ON (the `exceptions` preset, which also enables OSQP + NLopt)** — the exceptions
  carve-out: tests that assert the throwing wrappers throw (they need `REQUIRE_THROWS*`,
  unavailable under `CATCH_CONFIG_DISABLE_EXCEPTIONS`), the `[!shouldfail]` meta-test, and
  every OSQP/NLopt-linked test and bench, with Catch2 exceptions-on.

```sh
cmake --preset dev         # default: ctrlpp's tests + Catch2 under -fno-exceptions
cmake --preset exceptions  # carve-out: throwing-wrapper + OSQP/NLopt tests, exceptions-on
```

OSQP and NLopt throw internally and therefore require the exceptions build; the argmin static
NMPC path is throw-free and runs in the default `-fno-exceptions` tree.

## Documentation

- [Getting Started](docs/getting-started.md)<br/> Install ctrlpp and run your first PID controller
- [Guides](docs/guides/README.md)<br/> Tutorials and deep dives
- [API Reference](docs/README.md#api-reference)<br/> Full type documentation
- [Background Theory](docs/README.md#background-theory)<br/> Theory and mathematical background

## License

Apache 2.0 License, see [LICENSE](LICENSE) for details.

## Declaration of AI use
This library has been developed with extensive support from Claude Code in a hybrid of the Spiral Model and Extreme Programming using [GSD](https://github.com/gsd-build/get-shit-done).


