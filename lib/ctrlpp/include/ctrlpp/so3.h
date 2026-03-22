#ifndef HPP_GUARD_CTRLPP_SO3_H
#define HPP_GUARD_CTRLPP_SO3_H

// SO(3) Lie group primitives using unit quaternions.
//
// Convention: Hamilton convention throughout.
//   - Quaternion product q1*q2 corresponds to rotation q1 followed by q2.
//   - User-facing serialization is w-first: [w, x, y, z] (see to_vec/from_vec).
//   - Internally we use Eigen::Quaternion<Scalar>, which stores coefficients
//     in [x, y, z, w] order. We NEVER expose coeffs() directly to avoid the
//     w-last storage trap. Always use q.w(), q.vec(), q.x(), q.y(), q.z().
//
// References:
//   - Sola, "Quaternion kinematics for the error-state Kalman filter" (2017)
//   - Barfoot, "State Estimation for Robotics" (2017)

#include "ctrlpp/types.h"

#include <Eigen/Geometry>

#include <cmath>

namespace ctrlpp::so3 {

// Exponential map: rotation vector (angle-axis, phi) -> unit quaternion.
// Uses Rodrigues formula with Taylor expansion near zero to avoid division by zero.
template<typename Scalar>
[[nodiscard]] inline auto exp(const Vector<Scalar, 3>& phi) -> Eigen::Quaternion<Scalar>
{
    Scalar theta = phi.norm();
    Scalar half_theta = theta / Scalar{2};

    Scalar sinc_half;
    if (theta < Scalar{1e-7}) {
        // Taylor: sin(theta/2)/theta = 1/2 - theta^2/48 + ...
        sinc_half = Scalar{0.5} - theta * theta / Scalar{48};
    } else {
        sinc_half = std::sin(half_theta) / theta;
    }

    Eigen::Quaternion<Scalar> q;
    q.w() = std::cos(half_theta);
    q.vec() = sinc_half * phi;
    return q;
}

// Logarithmic map: unit quaternion -> rotation vector.
// Canonicalizes to w >= 0 hemisphere first for unique output.
template<typename Scalar>
[[nodiscard]] inline auto log(const Eigen::Quaternion<Scalar>& q) -> Vector<Scalar, 3>
{
    // Canonicalize: ensure w >= 0 (antipodal quaternions represent the same rotation)
    Eigen::Quaternion<Scalar> qc = q;
    if (qc.w() < Scalar{0}) {
        qc.w() = -qc.w();
        qc.x() = -qc.x();
        qc.y() = -qc.y();
        qc.z() = -qc.z();
    }

    Scalar vec_norm = qc.vec().norm();

    Scalar inv_sinc_half;
    if (vec_norm < Scalar{1e-7}) {
        // Taylor limit: 2 * atan2(eps, ~1) / eps -> 2
        inv_sinc_half = Scalar{2};
    } else {
        inv_sinc_half = Scalar{2} * std::atan2(vec_norm, qc.w()) / vec_norm;
    }

    return inv_sinc_half * qc.vec();
}

// Hamilton quaternion product: compose two rotations.
template<typename Scalar>
[[nodiscard]] inline auto compose(const Eigen::Quaternion<Scalar>& q1,
                                  const Eigen::Quaternion<Scalar>& q2) -> Eigen::Quaternion<Scalar>
{
    return q1 * q2;
}

// Quaternion conjugate (inverse for unit quaternions).
template<typename Scalar>
[[nodiscard]] inline auto conjugate(const Eigen::Quaternion<Scalar>& q) -> Eigen::Quaternion<Scalar>
{
    return q.conjugate();
}

// Normalize quaternion to unit norm.
template<typename Scalar>
[[nodiscard]] inline auto normalize(const Eigen::Quaternion<Scalar>& q) -> Eigen::Quaternion<Scalar>
{
    return q.normalized();
}

// Skew-symmetric matrix from a 3-vector: [v]_x such that [v]_x * u = v x u.
template<typename Scalar>
[[nodiscard]] inline auto skew(const Vector<Scalar, 3>& v) -> Matrix<Scalar, 3, 3>
{
    Matrix<Scalar, 3, 3> S;
    S <<  Scalar{0}, -v(2),  v(1),
           v(2),  Scalar{0}, -v(0),
          -v(1),   v(0), Scalar{0};
    return S;
}

// Quaternion to w-first vector: [w, x, y, z].
template<typename Scalar>
[[nodiscard]] inline auto to_vec(const Eigen::Quaternion<Scalar>& q) -> Vector<Scalar, 4>
{
    Vector<Scalar, 4> v;
    v << q.w(), q.vec();
    return v;
}

// W-first vector [w, x, y, z] to quaternion.
template<typename Scalar>
[[nodiscard]] inline auto from_vec(const Vector<Scalar, 4>& v) -> Eigen::Quaternion<Scalar>
{
    Eigen::Quaternion<Scalar> q;
    q.w() = v(0);
    q.vec() = v.template tail<3>();
    return q;
}

}

#endif
