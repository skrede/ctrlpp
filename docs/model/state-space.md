# state_space

Linear state-space model representations for continuous-time and discrete-time systems. The fundamental building block for control design, estimation, discretisation, and system identification throughout the library.

## Header and Alias

| Form | Header |
|------|--------|
| `continuous_state_space<Scalar, NX, NU, NY>` | `#include <ctrlpp/model/state_space.h>` |
| `discrete_state_space<Scalar, NX, NU, NY>` | `#include <ctrlpp/model/state_space.h>` |
| (convenience) | `#include <ctrlpp/state_space.h>` |

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | floating-point | Numeric type (`double`, `float`) |
| `NX` | `>= 1` | State dimension |
| `NU` | `>= 1` | Input dimension |
| `NY` | `>= 1` | Output dimension |

## continuous_state_space

Represents a continuous-time LTI system: `x_dot = A*x + B*u`, `y = C*x + D*u`.

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct continuous_state_space {
    Matrix<Scalar, NX, NX> A;
    Matrix<Scalar, NX, NU> B;
    Matrix<Scalar, NY, NX> C;
    Matrix<Scalar, NY, NU> D;
};
```

## discrete_state_space

Represents a discrete-time LTI system: `x[k+1] = A*x[k] + B*u[k]`, `y[k] = C*x[k] + D*u[k]`.

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct discrete_state_space {
    Matrix<Scalar, NX, NX> A;
    Matrix<Scalar, NX, NU> B;
    Matrix<Scalar, NY, NX> C;
    Matrix<Scalar, NY, NU> D;
};
```

## Type Aliases

```cpp
template <typename S, std::size_t NX>
using siso_continuous_state_space = continuous_state_space<S, NX, 1, 1>;

template <typename S, std::size_t NX>
using siso_discrete_state_space = discrete_state_space<S, NX, 1, 1>;
```

## Usage Example

```cpp
#include "ctrlpp/model/state_space.h"
#include "ctrlpp/model/analysis.h"
#include "ctrlpp/model/discretise.h"
#include "ctrlpp/model/propagate.h"

#include <Eigen/Dense>

#include <iostream>

int main()
{
    // Double integrator: position and velocity states, force input
    ctrlpp::continuous_state_space<double, 2, 1, 1> sys{
        .A = (Eigen::Matrix2d() << 0.0, 1.0, 0.0, 0.0).finished(),
        .B = (Eigen::Vector2d() << 0.0, 1.0).finished(),
        .C = (Eigen::RowVector2d() << 1.0, 0.0).finished(),
        .D = Eigen::Matrix<double, 1, 1>::Zero()};

    std::cout << "Continuous system:\n"
              << "  A =\n" << sys.A << "\n"
              << "  B = " << sys.B.transpose() << "\n"
              << "  Stable: " << ctrlpp::is_stable(sys) << "\n\n";

    // Discretise with ZOH at 100 Hz
    auto dsys = ctrlpp::discretise(sys, 0.01);

    // Simulate 10 steps with unit input
    Eigen::Vector2d x = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 1, 1> u;
    u << 1.0;

    for(int k = 0; k < 10; ++k)
    {
        auto y = ctrlpp::output(dsys, x, u);
        std::cout << "k=" << k << "  x=[" << x.transpose()
                  << "]  y=" << y[0] << "\n";
        x = ctrlpp::propagate(dsys, x, u);
    }
}
```

## See Also

- [transfer-function](transfer-function.md) -- transfer function representation
- [discretise](discretise.md) -- continuous-to-discrete conversion
- [conversion](conversion.md) -- TF to SS and SS to TF
- [analysis](analysis.md) -- stability, controllability, observability
- [propagate](propagate.md) -- state propagation utilities
