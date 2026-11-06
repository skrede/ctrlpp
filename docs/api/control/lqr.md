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
| `Scalar` | `ctrlpp_floating_scalar` | Floating-point type (e.g. `double`, `float`); the classes and `lqr_gain_continuous` state the constraint, and the Riccati-based functions inherit it from `dare` |
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
ctrlpp::expected<Eigen::Matrix<Scalar, int(NU), int(NX)>, dare_error>
lqr_gain(const Matrix<Scalar, NX, NX>& A,
         const Matrix<Scalar, NX, NU>& B,
         const Matrix<Scalar, NX, NX>& Q,
         const Matrix<Scalar, NU, NU>& R);
```

Computes the infinite-horizon LQR gain K = (R + B'PB)^{-1} B'PA where P is the stabilizing solution of the DARE.

**The gain returned is the solver's own -- `dare_result::K` -- not a recomputation from P.** The Riccati solve forms this gain inside the verification that accepts P: it applies the same dimension-times-unit-roundoff rank test to `R + B'PB` that the symplectic-matrix build applies to `R`, checks the result is finite, and forms it at the equilibrated scale. Re-forming it here would have neither guarantee. Without a rank threshold a rank-deficient `R + B'PB` yields a least-squares solution over the leading rank columns, returned as if it were the gain; and at the top of the range the caller-scale sum overflows on poses whose equilibrated answer is an ordinary normal number, where the rank-revealing solve returns a **zero** gain -- no feedback at all -- rather than refusing. Neither failure is visible to a caller who receives it. See [dare](dare.md) for the measured comparison.

**Rejections carry the Riccati solver's own `dare_error` verbatim.** The gain is a function of that solve and has no failure mode of its own, so it forwards the enumerator rather than restating the cause under a second name -- an unstabilizable pair, a singular state matrix and a non-converged factorization send the caller to fix three different things, and an empty result would have told them none of it. See [dare](dare.md) for the enumerators.

### lqr_gain (with cross-weight)

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
ctrlpp::expected<Eigen::Matrix<Scalar, int(NU), int(NX)>, dare_error>
lqr_gain(const Matrix<Scalar, NX, NX>& A,
         const Matrix<Scalar, NX, NU>& B,
         const Matrix<Scalar, NX, NX>& Q,
         const Matrix<Scalar, NU, NU>& R,
         const Matrix<Scalar, NX, NU>& N);
```

Infinite-horizon LQR gain with state-input cross-weight N: K = (R + B'PB)^{-1} (B'PA + N'). Forwards `dare_error` for the same reason, and returns the solver's own gain for the same reason. The cross-weight Riccati solve reduces the problem to standard form and then corrects the reduced problem's gain K' back to the posed problem's by `K = K' + R^{-1}N'`, using the `R^{-1}N'` it already formed for the reduction -- so the caller's gain costs one addition here and no second factorization.

### lqr_gain_continuous

```cpp
template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU,
          detail::care_solve_method Method = detail::sign_function_care_method>
ctrlpp::expected<Eigen::Matrix<Scalar, int(NU), int(NX)>, care_error>
lqr_gain_continuous(const Matrix<Scalar, NX, NX>& A,
                    const Matrix<Scalar, NX, NU>& B,
                    const Matrix<Scalar, NX, NX>& Q,
                    const Matrix<Scalar, NU, NU>& R,
                    Method method_tag = {});
```

Continuous-time gain `K = R^{-1} B' P` where P solves `A'P + PA - PBR^{-1}B'P + Q = 0`. Reports through `care_error`, forwarding the continuous solver's enumerator. `Method` is one of the three method tags [care](dare.md#care) takes (`detail::sign_function_care_method`, the default; `detail::schur_care_method`; `detail::balanced_schur_care_method`); the conditioning policy is always the default `detail::pivot_ratio_conditioning`.

**Both weightings are equilibrated first**, by the same common divisor `ctrlpp::care` applies: the smallest power of two at or above their largest entry. This surface does not go through `ctrlpp::care`; it forms `R^{-1}` and the Hamiltonian itself so it can reuse one factorization for the gain, which makes it a second continuous entry point. Without the same equilibration the two would answer different populations on the same arguments, and a caller cannot be expected to know which of them rescales.

Nothing is unscaled afterwards, and that is not an omission: the gain is homogeneous of degree **zero** in the common weight scale, since dividing `R` by `s` multiplies `R^{-1}` by `s` while dividing the solution by `s`. The gain formed from the equilibrated quantities is the caller's gain exactly. Only a solution would need multiplying back, and this surface does not return one. See the [numerical-behavior guide](../../guides/patterns/numerical-behavior.md) for what the equilibration buys: swept over eighteen decades of a common rescale in both directions, the equilibrated solver answers all 1,152 draws of both common-scale families, where the unequilibrated default answers 652 and 655.

It makes three rejections of its own before calling anything:

| Condition | Enumerator |
| --- | --- |
| a NaN or infinite `A`, `B`, `Q` or `R` | `care_error::non_finite_input` |
| a rank-deficient `R` | `care_error::singular_r` |
| a Hamiltonian that overflowed while being assembled | `care_error::non_finite_input` |

The middle one is worth stating explicitly, because this surface forms `R^{-1}` itself through an `LDLT` factorization rather than going through the Hamiltonian build, and **that factorization fails quietly**: its solve zeroes the rank-deficient directions instead of producing infinities. Without the rank test, a zero `R` would produce a finite `R^{-1}` of zeros, an entirely finite Hamiltonian describing a plant with no control authority, and a sign-function iteration that stagnates on it -- reported as `sign_function_stagnated`, which sends the caller to look at convergence rather than at the weighting they passed.

Past those three, every method makes one further refusal, and it applies to every solve rather than to malformed input:

| Condition | Enumerator | Methods |
| --- | --- | --- |
| the extracted solution does not satisfy the counted Riccati residual bound, or does not place the closed-loop spectrum strictly in the open left half-plane | `care_error::sign_function_stagnated` | the default sign-function tag |
| an acceptance magnitude cannot be resolved at the scalar type's range | `care_error::sign_function_stagnated` | the default sign-function tag |
| the same two conditions, on a matrix extracted by a Schur method | `care_error::unverified_solution` | `schur_care_method`, `balanced_schur_care_method` |

**That verification runs on every accepted solve, under every method tag.** It is not a fallback behind a cheaper check, and it is not a property of the default path: the solver reports success only for a matrix it has substituted back into the equation the caller posed. The consequence a caller sees is that a pose the solver cannot answer accurately is declined rather than answered, and the population that changes most is large common weight scales -- see [numerical behavior](../../guides/patterns/numerical-behavior.md) for which poses those are and why the boundary is where it is.

The two enumerators are separate because they say different things, and neither is a synonym for the other. `sign_function_stagnated` names a Newton iteration behavior: its documented cases are a non-finite scaling factor, a non-contracting step, an unresolvable magnitude, and a budget exhausted. A Schur method runs no Newton iteration, so reporting a failed Schur verification through it would send the caller to look at convergence when nothing converged or failed to. `unverified_solution` says only what happened: a candidate was extracted and could not be certified. A caller who meets it on one tag can try another.

See [dare](dare.md) for what `care_error::singular_u11` covers; the continuous enumerator has the identical shape as its discrete counterpart.

#### How far the enumerator can be relied on

**Branch on `has_value()`. Treat the enumerator as a diagnostic for a human reader rather than as a control signal.**

The tables above say which condition produces which enumerator, and they hold wherever the condition itself is decidable in the scalar type. They stop holding below the precision's resolution boundary -- the point at which the quantity deciding the problem falls under `sqrt(epsilon)` relative to the operands carrying it. Below it, *which* refusal the solve carries is decided by rounding rather than by the input.

Measured on one bit-identical input whose imaginary-axis mode is visible through `Q` at `4.2e-25` relative, a stabilizing solution that exists and is not determined in binary64, and on its thirty-point three-ulp neighborhood (`A(1,0)`, `A(1,1)`, `B(0)`, `B(1)` and `R` each moved by one to three units in the last place), over g++ 16.2.1 and clang 23.1.3, 21.1.8 and 18.1.8 on x86-64, `-O0` to `-O3`, fused multiply-add off and on, and Eigen 3.4.0 and 3.4.1:

| Input | Enumerator, over all 64 configurations |
| --- | --- |
| the input itself | `non_lhp_stabilizable` on every configuration |
| its thirty neighbors | all three of `non_lhp_stabilizable`, `non_psd_solution` and `sign_function_stagnated` on every configuration, zero acceptances; 27 of the 30 neighbors return different enumerators on different configurations |

arm64 and MSVC are not among the measured toolchains.

Sweeping the visibility across thirty-one decades, with the same neighborhood at each, puts the accept/refuse transition at `1e-8`, against `sqrt(epsilon) = 1.5e-8`: every neighbor is accepted down to `1e-7` on every configuration, and from `1e-8` to `1e-10` the verdict flips inside the neighborhood. `1e-11` and `1e-12` refuse all thirty neighbors on every configuration. `1e-13` accepts two of the thirty on 42 of the 64 configurations, every build without contraction or at `-O0` and g++ at `-O1` with it; under g++ at `-O2` those two answers agree with a binary128 Newton refinement to `1.3e-10` relative. Every decade from `1e-14` down refuses all thirty neighbors on every configuration while the enumerator keeps moving.

Two enumerators stay meaningful at any distance from the boundary, because they report a property of the operands rather than of an iteration: `non_finite_input` and `singular_r`. A finite, well-formed input with a nonsingular `R` never carries either.

#### Tuning the default tag

`sign_function_care_method` carries one defaulted member.

```cpp
struct sign_function_care_method
{
    sign_function_care_method(); // warmup_iterations{3}

    std::int32_t warmup_iterations;
};
```

**Neither `lqr_gain_continuous` nor `care` forwards the tag's value**: both take the tag for its type and run the solve with a default-constructed one, so `warmup_iterations` is `3` whatever the caller sets (see [Limitations](#limitations)). What follows describes the member as the iteration reads it.

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

Finite-horizon LQR via backward Riccati recursion from `P_N = Qf`. Returns gain sequence {K_0, K_1, ..., K_{N-1}} indexed by time step; `horizon = 0` gives an empty sequence. `Qf` is the terminal state cost. It has no failure channel (see [Limitations](#limitations)).

### lqr_tv_gains

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
ctrlpp::expected<std::vector<Eigen::Matrix<Scalar, int(NU), int(NX)>>, lqr_error>
lqr_tv_gains(const std::vector<Matrix<Scalar, NX, NX>>& As,
             const std::vector<Matrix<Scalar, NX, NU>>& Bs,
             const std::vector<Matrix<Scalar, NX, NX>>& Qs,
             const std::vector<Matrix<Scalar, NU, NU>>& Rs,
             const Matrix<Scalar, NX, NX>& Qf);
```

Time-varying LQR via backward Riccati recursion from `P_N = Qf` with per-step system and cost matrices. The horizon is the common length of the four sequences; empty sequences give an empty gain sequence. Refuses with:

- `lqr_error::sequence_length_mismatch` when the four sequences hold different counts.
- `lqr_error::non_finite_input` when any entry of any matrix, `Qf` included, is NaN or infinite.

### lqi_gain

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
ctrlpp::expected<lqi_result<Scalar, NX, NU, NY>, dare_error>
lqi_gain(const Matrix<Scalar, NX, NX>& A,
         const Matrix<Scalar, NX, NU>& B,
         const Matrix<Scalar, NY, NX>& C,
         const Matrix<Scalar, NX + NY, NX + NY>& Q_aug,
         const Matrix<Scalar, NU, NU>& R);
```

LQR with integral action. Augments the state with the discrete integral of the output, `A_aug = [[A, 0], [-C, I]]` and `B_aug = [[B], [0]]`, solves the augmented problem through `lqr_gain`, and returns `lqi_result` containing partitioned gains Kx (NU x NX) and Ki (NU x NY), or the augmented Riccati solve's `dare_error` verbatim. The augmentation carries no sample-time factor, so the integral weight in `Q_aug` absorbs it; `A_aug` is singular whenever `A` is, which `dare` refuses as `singular_a`.

### lqr_cost

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
ctrlpp::expected<Scalar, lqr_error>
lqr_cost(std::span<const Vector<Scalar, NX>> xs,
         std::span<const Vector<Scalar, NU>> us,
         const Matrix<Scalar, NX, NX>& Q,
         const Matrix<Scalar, NU, NU>& R);
```

Evaluates the quadratic trajectory cost: sum of x'Qx + u'Ru. `xs` holds either as many states as `us` holds inputs, or one more, in which case the terminal state cost is included. Refuses with:

- `lqr_error::trajectory_length_mismatch` for any other state count.
- `lqr_error::non_finite_input` when a state, an input, `Q` or `R` holds a NaN or an infinity.

## lqr Class

### Constructor

```cpp
explicit lqr(gain_type K);
```

Constructs a controller from a precomputed gain matrix.

### Methods

#### compute

```cpp
input_type compute(const state_type& x) const;
```

Returns u = -Kx.

#### gain

```cpp
const gain_type& gain() const;
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
ctrlpp::expected<input_type, lqr_error> compute(const state_type& x, std::size_t k) const;
```

Returns u = -K_k * x at time step k, or `lqr_error::step_out_of_range` when `k >= horizon()`.

#### gain

```cpp
ctrlpp::expected<gain_type, lqr_error> gain(std::size_t k) const;
```

Returns a copy of the gain matrix at step k, or `lqr_error::step_out_of_range` when `k >= horizon()`.

#### horizon

```cpp
std::size_t horizon() const;
```

Returns the number of time steps in the gain sequence.

## Supporting Types

### lqr_error

```cpp
enum class lqr_error
{
    sequence_length_mismatch,
    trajectory_length_mismatch,
    step_out_of_range,
    non_finite_input,
};
```

Declared in `ctrlpp/control/lqr_types.h` and included by `ctrlpp/control/lqr.h`. The refusals of `lqr_tv_gains`, `lqr_cost` and `lqr_time_varying::compute` and `gain`. Each is an exact precondition of the indexing or of the arithmetic's domain, so none carries a tolerance.

### lqi_result

```cpp
template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct lqi_result
{
    Eigen::Matrix<Scalar, int(NU), int(NX)> Kx;  // state feedback gain
    Eigen::Matrix<Scalar, int(NU), int(NY)> Ki;  // integral feedback gain
};
```

Declared in `ctrlpp/control/lqr_types.h` beside `lqr_error`, and included by `ctrlpp/control/lqr.h`.

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

## Limitations

### The finite-horizon recursions do not test their factorization

`lqr_finite` and `lqr_tv_gains` solve `(R + B'P B) K = B'P A` at every step with a column-pivoting QR factorization and no rank test, so a singular `R + B'P B` yields a solution of the rank-deficient system rather than a refusal. `lqr_finite` also tests nothing about its operands: a non-finite `A`, `B`, `Q`, `R` or `Qf` gives non-finite gains. `lqr_tv_gains` refuses non-finite operands and mismatched lengths but not a singular step.

### The sign-function tag's warmup_iterations is not forwarded

`lqr_gain_continuous` and both `care` overloads take the method tag as an argument but do not pass its value on: the sign-function solve always runs with `warmup_iterations = 3`. Setting the member on the tag the caller passes has no effect.

## See Also

- [dare](dare.md)<br/> discrete algebraic Riccati equation solver used internally
- [place](place.md)<br/> pole placement alternative to optimal control
- [kalman](../estimation/kalman.md)<br/> Kalman filter for observer-controller composition
- [guides/estimation/observer-controller](../../guides/estimation/observer-controller.md)<br/> observer-controller composition patterns
