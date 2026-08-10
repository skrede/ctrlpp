#ifndef HPP_GUARD_CTRLPP_LIE_SO3_H
#define HPP_GUARD_CTRLPP_LIE_SO3_H

/// @brief SO(3) Lie group primitives using unit quaternions.
///
/// Convention: Hamilton convention throughout.
///   - Quaternion product q1*q2 corresponds to rotation q1 followed by q2.
///   - User-facing serialization is w-first: [w, x, y, z] (see to_vec/from_vec).
///   - Internally we use Eigen::Quaternion<Scalar>, which stores coefficients
///     in [x, y, z, w] order. We NEVER expose coeffs() directly to avoid the
///     w-last storage trap. Always use q.w(), q.vec(), q.x(), q.y(), q.z().
///
/// @cite sola2018 -- Sola et al., "A micro Lie theory for state estimation in robotics", 2018
/// @cite barfoot2017 -- Barfoot, "State Estimation for Robotics", 2017, Ch. 7 (rotations, hemisphere canonicalisation of the log map)

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include <Eigen/Geometry>

#include <cmath>
#include <numbers>

namespace ctrlpp
{

/// @brief Structured failure modes for the SO(3) primitives that can refuse.
///
///  * non_finite_input : a quaternion coefficient is NaN or infinite, so no
///                       scaling of it lands on the unit sphere.
///  * zero_quaternion  : every coefficient is exactly zero. The quaternion
///                       carries no direction at all, and unlike a merely small
///                       one it cannot be brought onto the unit sphere by any
///                       factor.
///
/// The two are kept apart because they send the caller to fix different things:
/// the first is arithmetic that went wrong upstream, the second is a value that
/// was never a rotation.
enum class so3_error
{
    non_finite_input,
    zero_quaternion,
};

}

namespace ctrlpp::so3
{

/// Exponential map: rotation vector (angle-axis, phi) -> unit quaternion.
/// Uses Rodrigues formula with Taylor expansion near zero to avoid division by zero.
/// @cite sola2018 Eq. 101
template <ctrlpp_floating_scalar Scalar>
Eigen::Quaternion<Scalar> exp(const Vector<Scalar, 3>& phi)
{
    Scalar const scale = phi.cwiseAbs().maxCoeff();
    if(scale == Scalar{0})
        return Eigen::Quaternion<Scalar>::Identity();

    Vector<Scalar, 3> const scaled = phi / scale;
    Scalar const scaled_norm = scaled.norm();
    Scalar const theta = scale * scaled_norm;

    Eigen::Quaternion<Scalar> q;
    if(theta < static_cast<Scalar>(1e-7L))
    {
        Scalar const half_theta = theta / Scalar{2};
        Scalar const sinc_half = Scalar{0.5} - theta * theta / Scalar{48};
        q.w() = std::cos(half_theta);
        q.vec() = sinc_half * phi;
    }
    else
    {
        Scalar const period = Scalar{4} * std::numbers::pi_v<Scalar>;
        Scalar const scale_period = period / scaled_norm;
        Scalar const reduced_theta =
            std::fmod(std::fmod(scale, scale_period) * scaled_norm, period);
        Scalar const half_theta = reduced_theta / Scalar{2};
        q.w() = std::cos(half_theta);
        q.vec() = (std::sin(half_theta) / scaled_norm) * scaled;
    }

    return Eigen::Quaternion<Scalar>{q.coeffs() / q.norm()};
}

/// Logarithmic map: unit quaternion -> rotation vector.
/// Canonicalizes to w >= 0 hemisphere first for unique output.
/// @cite sola2018 Eq. 105
template <typename Scalar>
Vector<Scalar, 3> log(const Eigen::Quaternion<Scalar>& q)
{
    // Canonicalize: ensure w >= 0 (antipodal quaternions represent the same rotation)
    Eigen::Quaternion<Scalar> qc = q;
    if(qc.w() < Scalar{0})
    {
        qc.w() = -qc.w();
        qc.x() = -qc.x();
        qc.y() = -qc.y();
        qc.z() = -qc.z();
    }

    Scalar vec_norm = qc.vec().norm();

    Scalar inv_sinc_half;
    if(vec_norm < static_cast<Scalar>(1e-7))
        inv_sinc_half = Scalar{2}; // Taylor limit: 2 * atan2(eps, ~1) / eps -> 2
    else
        inv_sinc_half = Scalar{2} * std::atan2(vec_norm, qc.w()) / vec_norm;

    return inv_sinc_half * qc.vec();
}

/// Hamilton quaternion product: compose two rotations.
///
/// @cite sola2018 Sec. 4.3 (Hamilton quaternion product)
/// @cite barfoot2017 Ch. 7 (quaternion composition under Hamilton convention)
template <typename Scalar>
Eigen::Quaternion<Scalar> compose(const Eigen::Quaternion<Scalar>& q1, const Eigen::Quaternion<Scalar>& q2)
{
    return q1 * q2;
}

/// Quaternion conjugate (inverse for unit quaternions).
///
/// @cite sola2018 Sec. 4.3 (quaternion inverse via conjugation for unit quaternions)
template <typename Scalar>
Eigen::Quaternion<Scalar> conjugate(const Eigen::Quaternion<Scalar>& q)
{
    return q.conjugate();
}

/// Scale a quaternion onto the unit sphere, or report why it has no unit
/// representative.
///
/// The postcondition is a norm of one, and it holds for every input the
/// function accepts. Exactly two inputs have no unit representative and are
/// rejected rather than handed back: a NaN or infinite coefficient
/// (so3_error::non_finite_input), and the zero quaternion, which carries no
/// direction (so3_error::zero_quaternion).
///
/// Every other finite quaternion IS normalized, including ones whose squared
/// norm is not representable. The coefficients are divided by their largest
/// magnitude first, so the norm is taken of a vector whose largest coefficient
/// is exactly one and whose squared norm therefore lies in [1, 4] -- neither
/// end of the exponent range can be reached from there.
///
/// Forming the norm directly instead loses two families of finite input
/// silently, which is why this function does not delegate to the linear-algebra
/// library's normalizing member. That member tests squaredNorm() > 0 and
/// returns a COPY of its input when the test fails (Eigen 3.4.0,
/// Eigen/src/Core/Dot.h:122-134, MatrixBase<Derived>::normalized()), so a
/// quaternion whose squared norm underflows comes back unchanged with a norm of
/// zero; and where the squared norm overflows the test passes, the division is
/// by infinity, and the result is the zero quaternion. Both inputs are finite
/// and have a well-defined direction, and both would leave a documented
/// unit-norm postcondition unmet with nothing said about it.
template <typename Scalar>
auto normalize(const Eigen::Quaternion<Scalar>& q) -> ctrlpp::expected<Eigen::Quaternion<Scalar>, so3_error>
{
    if(!q.coeffs().allFinite())
        return ctrlpp::unexpected(so3_error::non_finite_input);

    Scalar const scale = q.coeffs().cwiseAbs().maxCoeff();
    if(!(scale > Scalar{0}))
        return ctrlpp::unexpected(so3_error::zero_quaternion);

    auto const scaled = (q.coeffs() / scale).eval();
    return Eigen::Quaternion<Scalar>{(scaled / scaled.norm()).eval()};
}

/// Skew-symmetric matrix from a 3-vector: [v]_x such that [v]_x * u = v x u.
/// @cite sola2018 Eq. 10
template <typename Scalar>
Matrix<Scalar, 3, 3> skew(const Vector<Scalar, 3>& v)
{
    Matrix<Scalar, 3, 3> S;
    S << Scalar{0}, -v(2), v(1), v(2), Scalar{0}, -v(0), -v(1), v(0), Scalar{0};
    return S;
}

// Quaternion to w-first vector: [w, x, y, z].
template <typename Scalar>
Vector<Scalar, 4> to_vec(const Eigen::Quaternion<Scalar>& q)
{
    Vector<Scalar, 4> v;
    v << q.w(), q.vec();
    return v;
}

// W-first vector [w, x, y, z] to quaternion.
template <typename Scalar>
Eigen::Quaternion<Scalar> from_vec(const Vector<Scalar, 4>& v)
{
    Eigen::Quaternion<Scalar> q;
    q.w() = v(0);
    q.vec() = v.template tail<3>();
    return q;
}

}

#endif
