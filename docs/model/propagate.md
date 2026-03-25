# propagate

State propagation and output computation utilities for discrete and continuous state-space models. Simple convenience functions that apply the state-space equations in a single call.

## Header and Alias

| Form | Header |
|------|--------|
| `propagate(sys, x, u)` | `#include <ctrlpp/model/propagate.h>` |
| `output(sys, x, u)` | `#include <ctrlpp/model/propagate.h>` |
| (convenience) | `#include <ctrlpp/propagate.h>` |

## Functions

### propagate

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
constexpr Vector<Scalar, NX>
propagate(const discrete_state_space<Scalar, NX, NU, NY>& sys,
          const Vector<Scalar, NX>& x,
          const Vector<Scalar, NU>& u);
```

Computes `x[k+1] = A*x[k] + B*u[k]`.

### output (discrete)

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
constexpr Vector<Scalar, NY>
output(const discrete_state_space<Scalar, NX, NU, NY>& sys,
       const Vector<Scalar, NX>& x,
       const Vector<Scalar, NU>& u);
```

Computes `y[k] = C*x[k] + D*u[k]`.

### output (continuous)

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
constexpr Vector<Scalar, NY>
output(const continuous_state_space<Scalar, NX, NU, NY>& sys,
       const Vector<Scalar, NX>& x,
       const Vector<Scalar, NU>& u);
```

Computes `y = C*x + D*u`.

## Usage Example

```cpp
// gnuplot: plot "< ./propagate_demo" using 1:2 with lines title "position"
#include <ctrlpp/model/propagate.h>
#include <ctrlpp/model/state_space.h>

#include <Eigen/Dense>

#include <iostream>

int main()
{
    constexpr double dt = 0.1;

    ctrlpp::discrete_state_space<double, 2, 1, 1> sys{
        .A = (Eigen::Matrix2d() << 1.0, dt, 0.0, 1.0).finished(),
        .B = (Eigen::Vector2d() << 0.5 * dt * dt, dt).finished(),
        .C = (Eigen::RowVector2d() << 1.0, 0.0).finished(),
        .D = Eigen::Matrix<double, 1, 1>::Zero()};

    Eigen::Vector2d x(0.0, 0.0);
    Eigen::Matrix<double, 1, 1> u;
    u << 1.0;

    for(int k = 0; k < 20; ++k)
    {
        auto y = ctrlpp::output(sys, x, u);
        std::cout << "k=" << k << "  pos=" << y[0]
                  << "  vel=" << x[1] << "\n";
        x = ctrlpp::propagate(sys, x, u);
    }
}
```

## See Also

- [state-space](state-space.md) -- state-space representation
- [dynamics-model](dynamics-model.md) -- nonlinear dynamics model concept
