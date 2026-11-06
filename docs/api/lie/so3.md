# so3

SO(3) Lie group primitives using Hamilton-convention unit quaternions. Provides free functions for exponential/logarithmic maps, quaternion composition, skew-symmetric matrix construction, and serialization. These primitives underpin the MEKF and manifold-UKF estimators.

Convention: Hamilton convention throughout, with `R(q1 * q2) = R(q1) R(q2)`: applied to a vector, `q1 * q2` rotates by `q2` and then by `q1` in the fixed frame, which read in the body frame is `q1` followed by `q2` about the axes `q1` produced. The attitude filters compose increments on the right, `q * exp(phi)`, so `phi` is a body-frame rotation vector. User-facing serialization is w-first: `[w, x, y, z]`. Internally, Eigen stores quaternion coefficients in `[x, y, z, w]` order.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::so3::*` (free functions), `ctrlpp::so3_error` | `#include <ctrlpp/lie/so3.h>` |
| (convenience) | `#include <ctrlpp/so3.h>` |

All functions live in the `ctrlpp::so3` namespace. Each takes a concrete `Vector<Scalar, 3>`, `Vector<Scalar, 4>` or `Eigen::Quaternion<Scalar>` and deduces `Scalar` from it, so an Eigen expression such as `omega * dt` is evaluated into a vector first (`Eigen::Vector3d step = omega * dt;` or `(omega * dt).eval()`); the expression type itself deduces nothing.

## Functions

### exp

```cpp
template <ctrlpp_floating_scalar Scalar>
Eigen::Quaternion<Scalar> exp(const Vector<Scalar, 3>& phi);
```

Exponential map: rotation vector (angle-axis) to unit quaternion, `[cos(theta/2), sin(theta/2) phi / theta]` with `theta = |phi|`. Below `theta = 1e-7` both parts take their Taylor series, `1 - theta^2/8` and `1/2 - theta^2/48`. The vector is divided by its largest component magnitude before its norm is formed, so a vector whose squared norm would overflow or underflow still maps correctly; the angle is reduced modulo `4 pi`, the period of the quaternion, and the result is normalized. The zero vector maps to the identity.

### log

```cpp
template <typename Scalar>
Vector<Scalar, 3> log(const Eigen::Quaternion<Scalar>& q);
```

Logarithmic map: unit quaternion to rotation vector, `2 atan2(|v|, w) v / |v|` for the vector part `v`. Because `q` and `-q` are the same rotation, the quaternion is first flipped to the `w >= 0` hemisphere, so `log(-q) == log(q)` and the norm, the rotation angle, lies in `[0, pi]`. At an angle of exactly `pi` (`w == 0`) nothing is flipped, so `q` and `-q` return the two opposite vectors of length `pi`, both of which name that rotation. Below `|v| = 1e-7` it returns `2 v`. The input is expected to be unit: the general branch is invariant to the quaternion's scale, the small-angle branch is not, and nothing is normalized or checked.

### compose

```cpp
template <typename Scalar>
Eigen::Quaternion<Scalar> compose(const Eigen::Quaternion<Scalar>& q1,
                                  const Eigen::Quaternion<Scalar>& q2);
```

Hamilton quaternion product, `q1 * q2`: composes two rotations under the convention above. The result is not renormalized.

### conjugate

```cpp
template <typename Scalar>
Eigen::Quaternion<Scalar> conjugate(const Eigen::Quaternion<Scalar>& q);
```

Quaternion conjugate (inverse for unit quaternions).

### normalize

```cpp
template <typename Scalar>
ctrlpp::expected<Eigen::Quaternion<Scalar>, so3_error> normalize(const Eigen::Quaternion<Scalar>& q);
```

Scales a quaternion onto the unit sphere, or reports why it has no unit
representative. The returned quaternion has a norm of one for every input the
function accepts.

Exactly two inputs are rejected, and they tell the caller to fix different
things:

| Enumerator | Condition | What it means |
| --- | --- | --- |
| `so3_error::non_finite_input` | a coefficient is NaN or infinite | arithmetic went wrong upstream; no scaling of this value lands on the unit sphere |
| `so3_error::zero_quaternion` | every coefficient is exactly zero | the value carries no direction and was never a rotation |

Every other finite quaternion is normalized, **including one whose squared norm
is not representable**. The coefficients are divided by their largest magnitude
before the norm is formed, so the norm is taken of a vector whose largest
coefficient is exactly one and whose squared norm lies in `[1, 4]`.

That is not a detail. Forming the norm directly loses two families of finite
input silently, which is why this function does not delegate to Eigen's
`normalized()`. That member tests `squaredNorm() > 0` and returns a **copy of
its input** when the test fails (Eigen 3.4.0,
`Eigen/src/Core/Dot.h:122-134`), so a quaternion whose squared norm underflows
comes back unchanged with a norm of zero; where the squared norm overflows
instead, the test passes, the division is by infinity, and the result is the
zero quaternion. Both inputs are finite and have a well-defined direction, so
both would leave a unit-norm postcondition unmet with nothing said about it.

```cpp
auto qn = ctrlpp::so3::normalize(q);
if(!qn.has_value())
{
    // qn.error() is so3_error::non_finite_input or so3_error::zero_quaternion
    return;
}
auto R = qn->toRotationMatrix();
```

### so3_error

```cpp
enum class so3_error
{
    non_finite_input,
    zero_quaternion,
};
```

Declared in namespace `ctrlpp` (not `ctrlpp::so3`), matching the other
per-module error enumerations.

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

Serializes a quaternion to a w-first 4-vector: `[w, x, y, z]`.

### from_vec

```cpp
template <typename Scalar>
Eigen::Quaternion<Scalar> from_vec(const Vector<Scalar, 4>& v);
```

Deserializes a w-first 4-vector `[w, x, y, z]` back to a quaternion. The coefficients are copied as given; the result is not normalized.

## Usage Example

```cpp
#include <ctrlpp/lie/so3.h>

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

- [estimation/mekf](../estimation/mekf.md)<br/> multiplicative EKF using SO(3)
- [estimation/manifold-ukf](../estimation/manifold-ukf.md)<br/> manifold UKF using SO(3)
- [background/so3](../../background/so3.md)<br/> SO(3), its Lie algebra and the exp/log maps
- [background/attitude-estimation](../../background/attitude-estimation.md)<br/> attitude estimation theory
