# discretize

Continuous-to-discrete state-space conversion. Converts a `continuous_state_space` to a `discrete_state_space` using zero-order hold (ZOH), Tustin (bilinear, with optional frequency prewarping), forward Euler, or backward Euler.

## Header and Alias

| Form | Header |
|------|--------|
| `discretize(sys, dt)` | `#include <ctrlpp/model/discretize.h>` |
| (convenience) | `#include <ctrlpp/discretize.h>` |

## Tag Types

```cpp
struct zoh {};            // Zero-order hold (default)
struct tustin {};         // Bilinear transform (Tustin)

template <typename Scalar>
struct tustin_prewarp     // Bilinear transform with frequency prewarping
{
    Scalar w_c;            // critical frequency, rad/s
};

struct forward_euler {};  // Forward Euler
struct backward_euler {}; // Backward Euler
```

## Functions

### discretize (ZOH, default)

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY>
discretize(const continuous_state_space<Scalar, NX, NU, NY>& sys,
           Scalar dt, zoh = {});
```

Discretizes using ZOH via the augmented matrix exponential. Forms the block matrix `[[A*dt, B*dt], [0, 0]]`, computes its matrix exponential, and extracts `Ad` and `Bd`. Output matrices `C` and `D` are passed through unchanged.

### discretize (explicit tag)

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY>
discretize(zoh, const continuous_state_space<Scalar, NX, NU, NY>& sys,
           Scalar dt);
```

Same as above with the tag as the first argument.

### discretize (Tustin / bilinear transform)

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY>
discretize(tustin, const continuous_state_space<Scalar, NX, NU, NY>& sys,
           Scalar dt);
```

Discretizes using the bilinear (Tustin) transform: `Ad = (I - A*dt/2)^-1 (I + A*dt/2)`, `Bd = (I - A*dt/2)^-1 B*dt`, `Cd = C (I - A*dt/2)^-1`. Includes the biproper feed-through correction `Dd = D + C*Bd/2`, which accounts for the direct coupling the bilinear map introduces between input and output even when the continuous system is strictly proper (`D = 0`).

### discretize (Tustin with frequency prewarping)

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY>
discretize(tustin_prewarp<Scalar> warp,
           const continuous_state_space<Scalar, NX, NU, NY>& sys, Scalar dt);
```

Rescales the sample period to `dt_warp = (2/w_c) * tan(w_c*dt/2)` before applying the bilinear map, so the discrete and continuous frequency responses agree exactly at the critical frequency `w_c` (rad/s). Everywhere else the pole mapping trades some accuracy for that exactness at `w_c`.

### discretize (forward Euler)

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY>
discretize(forward_euler, const continuous_state_space<Scalar, NX, NU, NY>& sys,
           Scalar dt);
```

Discretizes using the forward Euler approximation: `Ad = I + A*dt`, `Bd = B*dt`, `Cd = C`, `Dd = D`. First-order accurate; does not require a matrix inversion, but is only conditionally stable for a stable continuous system (the sample period must be small enough relative to the fastest pole).

### discretize (backward Euler)

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
discrete_state_space<Scalar, NX, NU, NY>
discretize(backward_euler, const continuous_state_space<Scalar, NX, NU, NY>& sys,
           Scalar dt);
```

Discretizes using the backward Euler approximation: `Ad = (I - A*dt)^-1`, `Bd = (I - A*dt)^-1 B*dt`, `Cd = C (I - A*dt)^-1`, `Dd = D + C (I - A*dt)^-1 B*dt`. First-order accurate and unconditionally stable for a stable continuous system, at the cost of a fixed-size matrix inversion.

## Usage Example

```cpp
// gnuplot: plot "< ./discretize_demo" using 1:2 with lines title "step response"
#include <ctrlpp/model/discretize.h>
#include <ctrlpp/model/state_space.h>
#include <ctrlpp/model/analysis.h>
#include <ctrlpp/model/propagate.h>

#include <Eigen/Dense>

#include <iostream>

int main()
{
    // Mass-spring-damper: m=1, k=4, b=1
    // x = [position, velocity], u = force
    ctrlpp::continuous_state_space<double, 2, 1, 1> sys{
        .A = (Eigen::Matrix2d() << 0.0, 1.0, -4.0, -1.0).finished(),
        .B = (Eigen::Vector2d() << 0.0, 1.0).finished(),
        .C = (Eigen::RowVector2d() << 1.0, 0.0).finished(),
        .D = Eigen::Matrix<double, 1, 1>::Zero()};

    std::cout << "Continuous stable: " << ctrlpp::is_stable(sys) << "\n";

    // Discretize with ZOH at 100 Hz
    constexpr double dt = 0.01;
    auto dsys = ctrlpp::discretize(sys, dt);

    std::cout << "Discrete A =\n" << dsys.A << "\n"
              << "Discrete B = " << dsys.B.transpose() << "\n"
              << "Discrete stable: " << ctrlpp::is_stable(dsys) << "\n\n";

    // Simulate step response
    Eigen::Vector2d x = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 1, 1> u;
    u << 1.0;

    for(int k = 0; k < 50; ++k)
    {
        auto y = ctrlpp::output(dsys, x, u);
        if(k % 10 == 0)
            std::cout << "t=" << k * dt << "  y=" << y[0] << "\n";
        x = ctrlpp::propagate(dsys, x, u);
    }
}
```

## See Also

- [state-space](state-space.md)<br/> state-space representations
- [propagate](propagate.md)<br/> propagate discretized systems
