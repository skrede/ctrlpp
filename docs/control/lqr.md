# lqr

Linear Quadratic Regulator providing infinite-horizon, finite-horizon, time-varying, and integral-action (LQI) gain computation, plus thin controller wrappers that store a precomputed gain matrix. The infinite-horizon variant solves the discrete algebraic Riccati equation (DARE) internally and returns the optimal state-feedback gain K such that u = -Kx minimises the quadratic cost J = sum(x'Qx + u'Ru).

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::lqr<Scalar, NX, NU>` | `#include <ctrlpp/control/lqr.h>` |
| `ctrlpp::lqr<Scalar, NX, NU>` | `#include <ctrlpp/lqr.h>` (convenience) |

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | arithmetic type | Numeric type (e.g. `double`, `float`) |
| `NX` | `std::size_t` | State dimension |
| `NU` | `std::size_t` | Input dimension |

## Type Aliases

```cpp
using gain_type  = Eigen::Matrix<Scalar, int(NU), int(NX)>;
using state_type = Eigen::Matrix<Scalar, int(NX), 1>;
using input_type = Eigen::Matrix<Scalar, int(NU), 1>;
```

## Free Functions

### lqr_gain

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
std::optional<Eigen::Matrix<Scalar, int(NU), int(NX)>>
lqr_gain(const Matrix<Scalar, NX, NX>& A,
         const Matrix<Scalar, NX, NU>& B,
         const Matrix<Scalar, NX, NX>& Q,
         const Matrix<Scalar, NU, NU>& R);
```

Computes the infinite-horizon LQR gain K = (R + B'PB)^{-1} B'PA where P is the stabilising solution of the DARE. Returns `std::nullopt` if the system is not stabilisable.

### lqr_gain (with cross-weight)

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
std::optional<Eigen::Matrix<Scalar, int(NU), int(NX)>>
lqr_gain(const Matrix<Scalar, NX, NX>& A,
         const Matrix<Scalar, NX, NU>& B,
         const Matrix<Scalar, NX, NX>& Q,
         const Matrix<Scalar, NU, NU>& R,
         const Matrix<Scalar, NX, NU>& N);
```

Infinite-horizon LQR gain with state-input cross-weight N: K = (R + B'PB)^{-1} (B'PA + N').

### lqr_finite

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
std::vector<Eigen::Matrix<Scalar, int(NU), int(NX)>>
lqr_finite(const Matrix<Scalar, NX, NX>& A,
           const Matrix<Scalar, NX, NU>& B,
           const Matrix<Scalar, NX, NX>& Q,
           const Matrix<Scalar, NU, NU>& R,
           const Matrix<Scalar, NX, NX>& Qf,
           std::size_t horizon);
```

Finite-horizon LQR via backward Riccati recursion. Returns gain sequence {K_0, K_1, ..., K_{N-1}} indexed by time step. `Qf` is the terminal state cost.

### lqr_tv_gains

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
std::vector<Eigen::Matrix<Scalar, int(NU), int(NX)>>
lqr_tv_gains(const std::vector<Matrix<Scalar, NX, NX>>& As,
             const std::vector<Matrix<Scalar, NX, NU>>& Bs,
             const std::vector<Matrix<Scalar, NX, NX>>& Qs,
             const std::vector<Matrix<Scalar, NU, NU>>& Rs,
             const Matrix<Scalar, NX, NX>& Qf,
             std::size_t horizon);
```

Time-varying LQR via backward Riccati recursion with per-step system and cost matrices.

### lqi_gain

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
std::optional<lqi_result<Scalar, NX, NU, NY>>
lqi_gain(const Matrix<Scalar, NX, NX>& A,
         const Matrix<Scalar, NX, NU>& B,
         const Matrix<Scalar, NY, NX>& C,
         const Matrix<Scalar, NX + NY, NX + NY>& Q_aug,
         const Matrix<Scalar, NU, NU>& R);
```

LQR with integral action. Augments the state with integral of tracking error and returns `lqi_result` containing partitioned gains Kx (NU x NX) and Ki (NU x NY).

### lqr_cost

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
Scalar lqr_cost(std::span<const Vector<Scalar, NX>> xs,
                std::span<const Vector<Scalar, NU>> us,
                const Matrix<Scalar, NX, NX>& Q,
                const Matrix<Scalar, NU, NU>& R);
```

Evaluates the quadratic trajectory cost: sum of x'Qx + u'Ru. If `xs` has one more element than `us`, the terminal state cost is included.

## lqr Class

### Constructor

```cpp
explicit lqr(gain_type K);
```

Constructs a controller from a precomputed gain matrix.

### Methods

#### compute

```cpp
auto compute(const state_type& x) const -> input_type;
```

Returns u = -Kx.

#### gain

```cpp
auto gain() const -> const gain_type&;
```

Returns a const reference to the stored gain matrix.

## lqr_time_varying Class

### Constructor

```cpp
explicit lqr_time_varying(std::vector<gain_type> gains);
```

Constructs a time-varying controller from a precomputed gain sequence.

### Methods

#### compute

```cpp
auto compute(const state_type& x, std::size_t k) const -> input_type;
```

Returns u = -K_k * x at time step k.

#### gain

```cpp
auto gain(std::size_t k) const -> const gain_type&;
```

Returns gain matrix at step k.

#### horizon

```cpp
auto horizon() const -> std::size_t;
```

Returns the number of time steps in the gain sequence.

## Supporting Types

### lqi_result

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct lqi_result
{
    Eigen::Matrix<Scalar, int(NU), int(NX)> Kx;  // state feedback gain
    Eigen::Matrix<Scalar, int(NU), int(NY)> Ki;  // integral feedback gain
};
```

## Usage Example

```cpp
#include <ctrlpp/control/lqr.h>
#include <ctrlpp/model/discretise.h>
#include <ctrlpp/model/state_space.h>

#include <Eigen/Dense>

#include <iostream>

int main()
{
    using Scalar = double;
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;

    // Mass-spring-damper: m=1, k=1, b=0.5, dt=0.05
    ctrlpp::continuous_state_space<Scalar, NX, NU, 1> sys_c{};
    sys_c.A << 0.0, 1.0, -1.0, -0.5;
    sys_c.B << 0.0, 1.0;
    sys_c.C << 1.0, 0.0;
    sys_c.D << 0.0;

    auto sys_d = ctrlpp::discretise(ctrlpp::zoh{}, sys_c, 0.05);

    Eigen::Matrix2d Q = Eigen::Matrix2d::Zero();
    Q(0, 0) = 10.0;
    Q(1, 1) = 1.0;
    Eigen::Matrix<Scalar, 1, 1> R;
    R << 1.0;

    auto K_opt = ctrlpp::lqr_gain<Scalar, NX, NU>(sys_d.A, sys_d.B, Q, R);
    if (!K_opt) {
        std::cerr << "DARE failed\n";
        return 1;
    }

    ctrlpp::lqr<Scalar, NX, NU> controller(*K_opt);

    Eigen::Vector2d x;
    x << 1.0, 0.0;

    for (int k = 0; k < 100; ++k) {
        auto u = controller.compute(x);
        x = sys_d.A * x + sys_d.B * u;
        std::cout << k * 0.05 << "," << x(0) << "," << u(0) << "\n";
    }
}
```

## See Also

- [dare](dare.md) -- discrete algebraic Riccati equation solver used internally
- [place](place.md) -- pole placement alternative to optimal control
- [kalman](../estimation/kalman.md) -- Kalman filter for observer-controller composition
- [guides/estimation/observer-controller](../guides/estimation/observer-controller.md) -- observer-controller composition patterns
