# analysis

System analysis functions for state-space models: poles, stability, controllability, observability, and closed-loop stability checks. Operates on both continuous and discrete systems.

## Header and Alias

| Form | Header |
|------|--------|
| `poles(sys)`, `is_stable(sys)`, etc. | `#include <ctrlpp/model/analysis.h>` |
| (convenience) | `#include <ctrlpp/analysis.h>` |

## Provability Semantics

Every predicate on this page answers whether the property is **provable** from the matrices it is given, not whether the property is merely "not disproved". An indeterminate input is therefore a negative answer: if any matrix a predicate reads carries a NaN or an infinity, nothing is provable from it and the predicate returns `false`. The [Limitations](#limitations) section states what a `false` does and does not tell the caller.

The guards are exact finiteness tests on the matrices each function reads. There is no tolerance and no threshold: a matrix either is entirely finite or it is not. `is_stable_closed_loop` and `is_stable_observer` additionally re-check the matrix they construct, because a finite operand triple can still overflow through the product, and an infinity minus an infinity leaves a NaN behind.

| Function | Matrices guarded | Answer for a non-finite input |
|----------|------------------|-------------------------------|
| `poles` | `sys.A` | all-NaN spectrum |
| `is_stable` (continuous) | `sys.A` | `false` |
| `is_stable` (discrete) | `sys.A` | `false` |
| `is_controllable` | `A`, `B` | `false` |
| `is_observable` | `A`, `C` | `false` |
| `is_stable_closed_loop` | `A`, `B`, `K`, and `A - B*K` | `false` |
| `is_stable_observer` | `A`, `L`, `C`, and `A - L*C` | `false` |

## Functions

### poles

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
std::array<std::complex<Scalar>, NX>
poles(const continuous_state_space<Scalar, NX, NU, NY>& sys);

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
std::array<std::complex<Scalar>, NX>
poles(const discrete_state_space<Scalar, NX, NU, NY>& sys);
```

Returns the eigenvalues of the A matrix (system poles).

The returned array is either every eigenvalue or explicitly none: a non-finite `sys.A`, and an eigen-solve that does not converge, both fill the whole array with a quiet NaN rather than leaving a partially populated one, so an unspecified solver output cannot be mistaken for a computed spectrum.

### is_stable

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
bool is_stable(const continuous_state_space<Scalar, NX, NU, NY>& sys);

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
bool is_stable(const discrete_state_space<Scalar, NX, NU, NY>& sys);
```

Continuous: `true` iff every computed pole has a negative real part. Discrete: `true` iff every computed pole has magnitude less than 1. Both comparisons are exact, with no margin. A non-finite `sys.A`, and a spectrum the solver could not compute, both answer `false`.

### is_controllable

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
bool is_controllable(const Matrix<Scalar, NX, NX>& A,
                     const Matrix<Scalar, NX, NU>& B);
```

Checks the rank of the controllability matrix `[B, AB, A^2 B, ..., A^{n-1} B]`, computed by `Eigen::FullPivLU` with its default threshold. Returns `true` if the rank equals NX. A non-finite entry in `A` or `B` returns `false`, because the rank of a matrix containing a NaN is not defined.

`NX` and `NU` are not deduced from the arguments: `Matrix` is an alias whose dimensions pass through `static_cast<int>`, so the call names them, as in `is_controllable<double, 2, 1>(A, B)`. The same holds for `is_observable`, `is_stable_closed_loop` and `is_stable_observer`.

### is_observable

```cpp
template <typename Scalar, std::size_t NX, std::size_t NY>
bool is_observable(const Matrix<Scalar, NX, NX>& A,
                   const Matrix<Scalar, NY, NX>& C);
```

Checks the rank of the observability matrix `[C; CA; CA^2; ...; CA^{n-1}]`, computed the same way. Returns `true` if the rank equals NX. A non-finite entry in `A` or `C` returns `false`, for the same reason as the controllability test.

### is_stable_closed_loop

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
bool is_stable_closed_loop(const Matrix<Scalar, NX, NX>& A,
                           const Matrix<Scalar, NX, NU>& B,
                           const Matrix<Scalar, NU, NX>& K);
```

Checks if all eigenvalues of `A - B*K` are inside the unit circle (discrete-time closed-loop stability; there is no continuous-time variant). Returns `false` when any of `A`, `B`, `K` is non-finite, when the constructed `A - B*K` is non-finite, which a finite operand triple can still produce through an overflowing product, and when the eigensolver does not converge.

### is_stable_observer

```cpp
template <typename Scalar, std::size_t NX, std::size_t NY>
bool is_stable_observer(const Matrix<Scalar, NX, NX>& A,
                        const Matrix<Scalar, NX, NY>& L,
                        const Matrix<Scalar, NY, NX>& C);
```

Checks if all eigenvalues of `A - L*C` are inside the unit circle (discrete-time observer stability; there is no continuous-time variant). Returns `false` when any of `A`, `L`, `C` is non-finite, when the constructed `A - L*C` is non-finite, and when the eigensolver does not converge.

## Usage Example

```cpp
#include <ctrlpp/model/analysis.h>
#include <ctrlpp/model/state_space.h>

#include <Eigen/Dense>

#include <complex>
#include <iostream>

int main()
{
    ctrlpp::continuous_state_space<double, 2, 1, 1> sys{
        (Eigen::Matrix2d() << 0.0, 1.0, -2.0, -3.0).finished(),
        (Eigen::Vector2d() << 0.0, 1.0).finished(),
        (Eigen::RowVector2d() << 1.0, 0.0).finished(),
        Eigen::Matrix<double, 1, 1>::Zero()};

    // Poles
    auto p = ctrlpp::poles(sys);
    std::cout << "Poles: ";
    for(const auto& pole : p)
        std::cout << pole << " ";
    std::cout << "\n";

    // Stability
    std::cout << "Stable: " << ctrlpp::is_stable(sys) << "\n";

    // Controllability and observability
    std::cout << "Controllable: "
              << ctrlpp::is_controllable<double, 2, 1>(sys.A, sys.B) << "\n";
    std::cout << "Observable: "
              << ctrlpp::is_observable<double, 2, 1>(sys.A, sys.C) << "\n";
}
```

## Limitations

### A false answer means not shown to hold

Every predicate returns a plain `bool`, with no third value for an input whose answer cannot be established. Such an input reads as `false`, so `false` means "not shown to hold", not "shown not to hold". The conditions that produce a `false` without deciding the property:

- any matrix the predicate reads has a NaN or infinite entry;
- `is_stable_closed_loop` and `is_stable_observer`: the constructed `A - B*K` or `A - L*C` is non-finite, or the eigensolver does not converge;
- `is_stable`: the spectrum is not computed (the eigensolver does not converge), so `poles` returns all NaN.

A caller that must tell an unstable, uncontrollable or unobservable system from an undecided one checks the matrices with `allFinite()` and inspects `poles` before reading the answer.

### Rank is decided within rounding

`is_controllable` and `is_observable` count the pivots of `Eigen::FullPivLU` whose magnitude exceeds `NX * eps` times the largest pivot, `eps` being the machine epsilon of `Scalar`. A system whose controllability or observability matrix is rank-deficient only by less than that threshold reads as `false`, although it is controllable or observable in exact arithmetic. With `A = diag(1, 2)` and `B = [b; 1]` in `double`, `b = 1e-14` answers `true` and `b = 1e-15` answers `false`. A `false` therefore also covers "too close to rank-deficient to resolve in `Scalar`".

### The stability tests have no margin

`is_stable` compares each computed pole exactly against the imaginary axis or the unit circle, and the closed-loop and observer tests compare exactly against the unit circle. A pole within rounding of the boundary can be reported either way: a continuous pole computed at `-1e-300` answers `true`. Neither test measures a stability margin.

## See Also

- [state-space](state-space.md)<br/> state-space representation
- [lqr](../control/lqr.md)<br/> LQR design requires controllability
- [place](../control/place.md)<br/> pole placement for control/observer design
