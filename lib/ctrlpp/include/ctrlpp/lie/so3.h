#ifndef HPP_GUARD_CTRLPP_LIE_SO3_H
#define HPP_GUARD_CTRLPP_LIE_SO3_H

// Hamilton convention throughout: q1 * q2 is rotation q1 followed by q2. Serialization is w-first,
// [w, x, y, z], while Eigen::Quaternion stores [x, y, z, w], so coeffs() is never exposed; the
// components are read through w(), vec(), x(), y() and z().
// @cite sola2018 -- Sola, Deray & Atchuthan, "A micro Lie theory for state estimation in robotics", arXiv:1812.01537v9, 2021, App. B, Eqs. 131 and 138, and App. B-A after Eq. 133

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include <Eigen/Geometry>

#include <cmath>
#include <numbers>

namespace ctrlpp
{

enum class so3_error
{
    non_finite_input,
    zero_quaternion,
};

}

namespace ctrlpp::so3
{

// Rodrigues formula, with its Taylor expansion near zero.
// @cite sola2018 -- Sola, Deray & Atchuthan, "A micro Lie theory for state estimation in robotics", arXiv:1812.01537v9, 2021, Eq. 132
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
        Scalar const sinc_half = Scalar{0.5} - theta * theta / Scalar{48};
        q.w() = Scalar{1} - theta * theta / Scalar{8};
        q.vec() = sinc_half * phi;
    }
    else
    {
        Scalar const period = Scalar{4} * std::numbers::pi_v<Scalar>;
        Scalar const scale_period = period / scaled_norm;
        Scalar const reduced_theta = std::fmod(std::fmod(scale, scale_period) * scaled_norm, period);
        Scalar const half_theta = reduced_theta / Scalar{2};
        q.w() = std::cos(half_theta);
        q.vec() = (std::sin(half_theta) / scaled_norm) * scaled;
    }

    return Eigen::Quaternion<Scalar>{q.coeffs() / q.norm()};
}

// Canonicalizes to the w >= 0 hemisphere first, since antipodal quaternions are the same rotation
// and the output must be unique.
// @cite sola2018 -- Sola, Deray & Atchuthan, "A micro Lie theory for state estimation in robotics", arXiv:1812.01537v9, 2021, Eq. 133
template <typename Scalar>
Vector<Scalar, 3> log(const Eigen::Quaternion<Scalar>& q)
{
    Eigen::Quaternion<Scalar> qc = q;
    if(qc.w() < Scalar{0})
    {
        qc.w() = -qc.w();
        qc.x() = -qc.x();
        qc.y() = -qc.y();
        qc.z() = -qc.z();
    }

    Scalar vec_norm = qc.vec().norm();

    // Below the threshold, the Taylor limit of 2 * atan2(v, w) / v as v -> 0.
    Scalar inv_sinc_half;
    if(vec_norm < static_cast<Scalar>(1e-7))
        inv_sinc_half = Scalar{2};
    else
        inv_sinc_half = Scalar{2} * std::atan2(vec_norm, qc.w()) / vec_norm;

    return inv_sinc_half * qc.vec();
}

// @cite sola2018 -- Sola, Deray & Atchuthan, "A micro Lie theory for state estimation in robotics", arXiv:1812.01537v9, 2021, App. B, Eq. 131
template <typename Scalar>
Eigen::Quaternion<Scalar> compose(const Eigen::Quaternion<Scalar>& q1, const Eigen::Quaternion<Scalar>& q2)
{
    return q1 * q2;
}

// @cite sola2018 -- Sola, Deray & Atchuthan, "A micro Lie theory for state estimation in robotics", arXiv:1812.01537v9, 2021, App. B
template <typename Scalar>
Eigen::Quaternion<Scalar> conjugate(const Eigen::Quaternion<Scalar>& q)
{
    return q.conjugate();
}

// Divides by the largest coefficient magnitude before forming the norm, so a squared norm that would
// underflow or overflow still normalizes. Eigen's normalized() returns its input unchanged when
// squaredNorm() underflows and the zero quaternion when it overflows (Eigen 3.4.0,
// Eigen/src/Core/Dot.h:122-134, MatrixBase<Derived>::normalized()).
template <typename Scalar>
ctrlpp::expected<Eigen::Quaternion<Scalar>, so3_error> normalize(const Eigen::Quaternion<Scalar>& q)
{
    if(!q.coeffs().allFinite())
        return ctrlpp::unexpected(so3_error::non_finite_input);

    Scalar const scale = q.coeffs().cwiseAbs().maxCoeff();
    if(!(scale > Scalar{0}))
        return ctrlpp::unexpected(so3_error::zero_quaternion);

    auto const scaled = (q.coeffs() / scale).eval();
    return Eigen::Quaternion<Scalar>{(scaled / scaled.norm()).eval()};
}

// [v]_x * u = v x u.
// @cite sola2018 Eq. 10
template <typename Scalar>
Matrix<Scalar, 3, 3> skew(const Vector<Scalar, 3>& v)
{
    Matrix<Scalar, 3, 3> S;
    S << Scalar{0}, -v(2), v(1), v(2), Scalar{0}, -v(0), -v(1), v(0), Scalar{0};
    return S;
}

template <typename Scalar>
Vector<Scalar, 4> to_vec(const Eigen::Quaternion<Scalar>& q)
{
    Vector<Scalar, 4> v;
    v << q.w(), q.vec();
    return v;
}

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
