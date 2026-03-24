# analysis

System analysis functions for state-space models: poles, stability, controllability, observability, and closed-loop stability checks. Operates on both continuous and discrete systems.

## Header and Alias

| Form | Header |
|------|--------|
| `poles(sys)`, `is_stable(sys)`, etc. | `#include <ctrlpp/model/analysis.h>` |
| (convenience) | `#include <ctrlpp/analysis.h>` |

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

### is_stable

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
bool is_stable(const continuous_state_space<Scalar, NX, NU, NY>& sys);

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
bool is_stable(const discrete_state_space<Scalar, NX, NU, NY>& sys);
```

Continuous: stable iff all poles have negative real part. Discrete: stable iff all poles have magnitude less than 1.

### is_controllable

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
bool is_controllable(const Matrix<Scalar, NX, NX>& A,
                     const Matrix<Scalar, NX, NU>& B);
```

Checks rank of the controllability matrix `[B, AB, A^2 B, ..., A^{n-1} B]`. Returns `true` if rank equals NX.

### is_observable

```cpp
template <typename Scalar, std::size_t NX, std::size_t NY>
bool is_observable(const Matrix<Scalar, NX, NX>& A,
                   const Matrix<Scalar, NY, NX>& C);
```

Checks rank of the observability matrix `[C; CA; CA^2; ...; CA^{n-1}]`. Returns `true` if rank equals NX.

### is_stable_closed_loop

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
bool is_stable_closed_loop(const Matrix<Scalar, NX, NX>& A,
                           const Matrix<Scalar, NX, NU>& B,
                           const Matrix<Scalar, NU, NX>& K);
```

Checks if all eigenvalues of `A - B*K` are inside the unit circle (discrete-time closed-loop stability).

### is_stable_observer

```cpp
template <typename Scalar, std::size_t NX, std::size_t NY>
bool is_stable_observer(const Matrix<Scalar, NX, NX>& A,
                        const Matrix<Scalar, NX, NY>& L,
                        const Matrix<Scalar, NY, NX>& C);
```

Checks if all eigenvalues of `A - L*C` are inside the unit circle (discrete-time observer stability).

## Usage Example

```cpp
#include "ctrlpp/model/analysis.h"
#include "ctrlpp/model/state_space.h"

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

- [state-space](state-space.md) -- state-space representation
