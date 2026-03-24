# so3

SO(3) Lie group primitives using Hamilton-convention unit quaternions. Provides free functions for exponential/logarithmic maps, quaternion composition, skew-symmetric matrix construction, and serialisation. These primitives underpin the MEKF and manifold-UKF estimators.

Convention: Hamilton convention throughout. Quaternion product `q1 * q2` corresponds to rotation `q1` followed by `q2`. User-facing serialisation is w-first: `[w, x, y, z]`. Internally, Eigen stores quaternion coefficients in `[x, y, z, w]` order.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::so3::*` (free functions) | `#include <ctrlpp/lie/so3.h>` |
| (convenience) | `#include <ctrlpp/so3.h>` |

All functions live in the `ctrlpp::so3` namespace.

## Functions

### exp

```cpp
template <typename Scalar>
Eigen::Quaternion<Scalar> exp(const Vector<Scalar, 3>& phi);
```

Exponential map: rotation vector (angle-axis) to unit quaternion. Uses Rodrigues formula with Taylor expansion near zero for numerical stability.

### log

```cpp
template <typename Scalar>
Vector<Scalar, 3> log(const Eigen::Quaternion<Scalar>& q);
```

Logarithmic map: unit quaternion to rotation vector. Canonicalises to the `w >= 0` hemisphere for unique output.

### compose

```cpp
template <typename Scalar>
Eigen::Quaternion<Scalar> compose(const Eigen::Quaternion<Scalar>& q1,
                                  const Eigen::Quaternion<Scalar>& q2);
```

Hamilton quaternion product: composes two rotations.

### conjugate

```cpp
template <typename Scalar>
Eigen::Quaternion<Scalar> conjugate(const Eigen::Quaternion<Scalar>& q);
```

Quaternion conjugate (inverse for unit quaternions).

### normalize

```cpp
template <typename Scalar>
Eigen::Quaternion<Scalar> normalize(const Eigen::Quaternion<Scalar>& q);
```

Normalises a quaternion to unit norm.

### skew

```cpp
template <typename Scalar>
Matrix<Scalar, 3, 3> skew(const Vector<Scalar, 3>& v);
```

Constructs the skew-symmetric matrix `[v]_x` such that `[v]_x * u = v x u` (cross product).

### to_vec

```cpp
template <typename Scalar>
Vector<Scalar, 4> to_vec(const Eigen::Quaternion<Scalar>& q);
```

Serialises a quaternion to a w-first 4-vector: `[w, x, y, z]`.

### from_vec

```cpp
template <typename Scalar>
Eigen::Quaternion<Scalar> from_vec(const Vector<Scalar, 4>& v);
```

Deserialises a w-first 4-vector `[w, x, y, z]` back to a quaternion.

## Usage Example

```cpp
#include "ctrlpp/lie/so3.h"

#include <Eigen/Geometry>

#include <cmath>
#include <iostream>
#include <numbers>

int main()
{
    using namespace ctrlpp;

    // Create a 90-degree rotation about the z-axis
    Vector<double, 3> phi{0.0, 0.0, std::numbers::pi / 2.0};
    auto q1 = so3::exp(phi);

    std::cout << "q1 (90 deg about z): " << so3::to_vec(q1).transpose() << "\n";

    // Verify exp/log round-trip
    auto phi_back = so3::log(q1);
    std::cout << "log(q1): " << phi_back.transpose() << "\n";
    std::cout << "Round-trip error: " << (phi - phi_back).norm() << "\n\n";

    // Compose two rotations: 45 deg + 45 deg about z = 90 deg
    Vector<double, 3> half_phi{0.0, 0.0, std::numbers::pi / 4.0};
    auto q_half = so3::exp(half_phi);
    auto q_composed = so3::compose(q_half, q_half);

    std::cout << "45+45 composed: " << so3::to_vec(q_composed).transpose() << "\n";
    std::cout << "Difference from 90: "
              << so3::log(so3::compose(so3::conjugate(q1), q_composed)).norm() << "\n\n";

    // Skew-symmetric matrix
    Vector<double, 3> v{1.0, 2.0, 3.0};
    auto S = so3::skew(v);
    std::cout << "skew([1,2,3]):\n" << S << "\n";
}
```

## See Also

- [estimation/mekf](../estimation/mekf.md) -- multiplicative EKF using SO(3)
- [estimation/manifold-ukf](../estimation/manifold-ukf.md) -- manifold UKF using SO(3)
- [reference/attitude-estimation-theory](../reference/attitude-estimation-theory.md) -- attitude estimation theory
