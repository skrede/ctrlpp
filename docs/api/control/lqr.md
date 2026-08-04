# lqr

Linear Quadratic Regulator providing infinite-horizon, finite-horizon, time-varying, and integral-action (LQI) gain computation, plus thin controller wrappers that store a precomputed gain matrix. The infinite-horizon variant solves the discrete algebraic Riccati equation (DARE) internally and returns the optimal state-feedback gain K such that u = -Kx minimizes the quadratic cost J = sum(x'Qx + u'Ru).

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
auto lqr_gain(const Matrix<Scalar, NX, NX>& A,
              const Matrix<Scalar, NX, NU>& B,
              const Matrix<Scalar, NX, NX>& Q,
              const Matrix<Scalar, NU, NU>& R)
    -> ctrlpp::expected<Eigen::Matrix<Scalar, int(NU), int(NX)>, dare_error>;
```

Computes the infinite-horizon LQR gain K = (R + B'PB)^{-1} B'PA where P is the stabilizing solution of the DARE.

**Rejections carry the Riccati solver's own `dare_error` verbatim.** The gain is a function of that solve and has no failure mode of its own, so it forwards the enumerator rather than restating the cause under a second name -- an unstabilizable pair, a singular state matrix and a non-converged factorization send the caller to fix three different things, and an empty result would have told them none of it. See [dare](dare.md) for the enumerators.

### lqr_gain (with cross-weight)

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
auto lqr_gain(const Matrix<Scalar, NX, NX>& A,
              const Matrix<Scalar, NX, NU>& B,
              const Matrix<Scalar, NX, NX>& Q,
              const Matrix<Scalar, NU, NU>& R,
              const Matrix<Scalar, NX, NU>& N)
    -> ctrlpp::expected<Eigen::Matrix<Scalar, int(NU), int(NX)>, dare_error>;
```

Infinite-horizon LQR gain with state-input cross-weight N: K = (R + B'PB)^{-1} (B'PA + N'). Forwards `dare_error` for the same reason.

### lqr_gain_continuous

```cpp
template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU,
          detail::care_solve_method Method = detail::sign_function_care_method>
auto lqr_gain_continuous(const Matrix<Scalar, NX, NX>& A,
                         const Matrix<Scalar, NX, NU>& B,
                         const Matrix<Scalar, NX, NX>& Q,
                         const Matrix<Scalar, NU, NU>& R,
                         Method method_tag = {})
    -> ctrlpp::expected<Eigen::Matrix<Scalar, int(NU), int(NX)>, care_error>;
```

Continuous-time gain `K = R^{-1} B' P` where P solves `A'P + PA - PBR^{-1}B'P + Q = 0`. Reports through `care_error`, forwarding the continuous solver's enumerator.

It makes three rejections of its own before calling anything:

| Condition | Enumerator |
| --- | --- |
| a NaN or infinite `A`, `B`, `Q` or `R` | `care_error::non_finite_input` |
| a rank-deficient `R` | `care_error::singular_r` |
| a Hamiltonian that overflowed while being assembled | `care_error::non_finite_input` |

The middle one is worth stating explicitly, because this surface forms `R^{-1}` itself through an `LDLT` factorization rather than going through the Hamiltonian build, and **that factorization fails quietly**: its solve zeroes the rank-deficient directions instead of producing infinities. Before the rank test, a zero `R` therefore produced a finite `R^{-1}` of zeros, an entirely finite Hamiltonian describing a plant with no control authority, and a sign-function iteration that stagnated on it -- reported as `sign_function_stagnated`, which sends the caller to look at convergence rather than at the weighting they passed.

Past those three, every method makes one further refusal, and it applies to every solve rather than to malformed input:

| Condition | Enumerator | Methods |
| --- | --- | --- |
| the extracted solution does not satisfy the counted Riccati residual bound, or does not place the closed-loop spectrum strictly in the open left half-plane | `care_error::sign_function_stagnated` | the default sign-function tag |
| an acceptance magnitude cannot be resolved at the scalar type's range | `care_error::sign_function_stagnated` | the default sign-function tag |
| the same two conditions, on a matrix extracted by a Schur method | `care_error::unverified_solution` | `schur_care_method`, `balanced_schur_care_method` |

**That verification runs on every accepted solve, under every method tag.** It is not a fallback behind a cheaper check, and it is not a property of the default path: the solver reports success only for a matrix it has substituted back into the equation the caller posed. The consequence a caller sees is that a pose the solver cannot answer accurately is declined rather than answered, and the population that changes most is large common weight scales -- see [numerical behavior](../../guides/patterns/numerical-behavior.md) for which poses those are and why the boundary is where it is.

The two enumerators are separate because they say different things, and neither is a synonym for the other. `sign_function_stagnated` names a Newton iteration behavior: its documented cases are a non-finite scaling factor, a non-contracting step, an unresolvable magnitude, and a budget exhausted. A Schur method runs no Newton iteration, so reporting a failed Schur verification through it would send the caller to look at convergence when nothing converged or failed to. `unverified_solution` says only what happened: a candidate was extracted and could not be certified. A caller who meets it on one tag should try another -- the balanced variant in particular answers common weight rescales that the other two decline.

See [dare](dare.md) for what `care_error::singular_u11` covers; the continuous enumerator has the identical shape as its discrete counterpart.

#### Tuning the default tag

`sign_function_care_method` carries one defaulted member. Omitting the tag, or passing a default-constructed one, gives today's behavior exactly.

```cpp
struct sign_function_care_method
{
    int warmup_iterations = 3;
};

// Give a badly scaled Hamiltonian more room before the guard arms.
auto K = ctrlpp::lqr_gain_continuous<double, 2, 2>(
    A, B, Q, R, ctrlpp::detail::sign_function_care_method{.warmup_iterations = 8});
```

`warmup_iterations` is the number of Newton steps taken before the non-contraction guard arms. After the window, a step whose change grows rather than shrinks ends the solve with `care_error::sign_function_stagnated`; inside it, growth is allowed, because the determinantal scaling makes large corrections in the early steps. Raising it trades a later decline for a chance at an answer; lowering it declines sooner.

It is a knob rather than a constant because the number of such early steps is a property of the input's conditioning and is not derivable from the scalar type or the dimension. The default's provenance is stated rather than implied: swept over 3,456 draws spanning eighteen decades of weight scale in each direction, every value from 0 to the iteration cap of 40 produced identical outcomes, so `3` is retained because it moves nothing on that evidence.

**The knob cannot break the result contract.** The guard bounds wasted work; it does not decide acceptance. Whatever the iteration produces is still verified against the caller's own Hamiltonian before it is reported, under every value.

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
auto lqi_gain(const Matrix<Scalar, NX, NX>& A,
              const Matrix<Scalar, NX, NU>& B,
              const Matrix<Scalar, NY, NX>& C,
              const Matrix<Scalar, NX + NY, NX + NY>& Q_aug,
              const Matrix<Scalar, NU, NU>& R)
    -> ctrlpp::expected<lqi_result<Scalar, NX, NU, NY>, dare_error>;
```

LQR with integral action. Augments the state with integral of tracking error and returns `lqi_result` containing partitioned gains Kx (NU x NX) and Ki (NU x NY), or the augmented Riccati solve's `dare_error` verbatim.

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
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'position', '' using 1:3 with lines title 'control'"

#include <ctrlpp/control/lqr.h>
#include <ctrlpp/model/discretize.h>
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

    auto sys_d = ctrlpp::discretize(ctrlpp::zoh{}, sys_c, 0.05);

    Eigen::Matrix2d Q = Eigen::Matrix2d::Zero();
    Q(0, 0) = 10.0;
    Q(1, 1) = 1.0;
    Eigen::Matrix<Scalar, 1, 1> R;
    R << 1.0;

    auto K_opt = ctrlpp::lqr_gain<Scalar, NX, NU>(sys_d.A, sys_d.B, Q, R);
    if (!K_opt.has_value()) {
        // K_opt.error() is the dare_error naming which condition refused the
        // plant, not merely that something did.
        std::cerr << "the Riccati solve refused the plant\n";
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

- [dare](dare.md)<br/> discrete algebraic Riccati equation solver used internally
- [place](place.md)<br/> pole placement alternative to optimal control
- [kalman](../estimation/kalman.md)<br/> Kalman filter for observer-controller composition
- [guides/estimation/observer-controller](../../guides/estimation/observer-controller.md)<br/> observer-controller composition patterns
