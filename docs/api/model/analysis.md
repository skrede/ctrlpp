# analysis

System analysis functions for state-space models: poles, stability, controllability, observability, and closed-loop stability checks. Operates on both continuous and discrete systems.

## Header and Alias

| Form | Header |
|------|--------|
| `poles(sys)`, `is_stable(sys)`, etc. | `#include <ctrlpp/model/analysis.h>` |
| (convenience) | `#include <ctrlpp/analysis.h>` |

## Provability Semantics

Every predicate on this page answers whether the property is **provable** from the matrices it is given, not whether the property is merely "not disproved". An indeterminate input is therefore a negative answer: if any matrix a predicate reads carries a NaN or an infinity, nothing is provable from it and the predicate returns `false`.

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

Continuous: provably stable iff all poles have negative real part. Discrete: provably stable iff all poles have magnitude less than 1. A non-finite `sys.A`, and a spectrum the solver could not compute, both answer `false`.

### is_controllable

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
bool is_controllable(const Matrix<Scalar, NX, NX>& A,
                     const Matrix<Scalar, NX, NU>& B);
```

Checks rank of the controllability matrix `[B, AB, A^2 B, ..., A^{n-1} B]`. Returns `true` if rank equals NX. A non-finite entry in `A` or `B` returns `false`, because the rank of a matrix containing a NaN is not defined.

### is_observable

```cpp
template <typename Scalar, std::size_t NX, std::size_t NY>
bool is_observable(const Matrix<Scalar, NX, NX>& A,
                   const Matrix<Scalar, NY, NX>& C);
```

Checks rank of the observability matrix `[C; CA; CA^2; ...; CA^{n-1}]`. Returns `true` if rank equals NX. A non-finite entry in `A` or `C` returns `false`, for the same reason as the controllability test.

### is_stable_closed_loop

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
bool is_stable_closed_loop(const Matrix<Scalar, NX, NX>& A,
                           const Matrix<Scalar, NX, NU>& B,
                           const Matrix<Scalar, NU, NX>& K);
```

Checks if all eigenvalues of `A - B*K` are inside the unit circle (discrete-time closed-loop stability). Returns `false` when any of `A`, `B`, `K` is non-finite, and also when the constructed `A - B*K` is non-finite, which a finite operand triple can still produce through an overflowing product.

### is_stable_observer

```cpp
template <typename Scalar, std::size_t NX, std::size_t NY>
bool is_stable_observer(const Matrix<Scalar, NX, NX>& A,
                        const Matrix<Scalar, NX, NY>& L,
                        const Matrix<Scalar, NY, NX>& C);
```

Checks if all eigenvalues of `A - L*C` are inside the unit circle (discrete-time observer stability). Returns `false` when any of `A`, `L`, `C` is non-finite, and also when the constructed `A - L*C` is non-finite.

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
        .A = (Eigen::Matrix2d() << 0.0, 1.0, -2.0, -3.0).finished(),
        .B = (Eigen::Vector2d() << 0.0, 1.0).finished(),
        .C = (Eigen::RowVector2d() << 1.0, 0.0).finished(),
        .D = Eigen::Matrix<double, 1, 1>::Zero()};

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
              << ctrlpp::is_controllable(sys.A, sys.B) << "\n";
    std::cout << "Observable: "
              << ctrlpp::is_observable(sys.A, sys.C) << "\n";
}
```

## See Also

- [state-space](state-space.md)<br/> state-space representation
- [lqr](../control/lqr.md)<br/> LQR design requires controllability
- [place](../control/place.md)<br/> pole placement for control/observer design
