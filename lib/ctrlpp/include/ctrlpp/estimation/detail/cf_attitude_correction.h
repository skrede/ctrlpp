#ifndef HPP_GUARD_CTRLPP_ESTIMATION_DETAIL_CF_ATTITUDE_CORRECTION_H
#define HPP_GUARD_CTRLPP_ESTIMATION_DETAIL_CF_ATTITUDE_CORRECTION_H

// @cite mahony2008 -- Mahony et al., "Nonlinear Complementary Filters on the Special Orthogonal Group", 2008

#include "ctrlpp/types.h"

#include <Eigen/Geometry>

#include <cmath>

namespace ctrlpp::detail
{

// Scaled by the largest magnitude before the norm, so every finite nonzero reading has a direction even
// where its squared norm would underflow or overflow. No magnitude floor: only the scale-free direction is used.
template <typename Scalar>
inline bool unit_direction(const Vector<Scalar, 3>& v, Vector<Scalar, 3>& direction)
{
    const Scalar scale = v.cwiseAbs().maxCoeff();
    if(!(scale > Scalar{0}))
        return false;
    const Vector<Scalar, 3> scaled = (v / scale).eval();
    direction = (scaled / scaled.norm()).eval();
    return true;
}

// @cite mahony2008 -- Mahony et al., 2008, Eq. 12 (gravity error via cross product)
template <typename Scalar>
inline Vector<Scalar, 3> compute_gravity_correction(const Eigen::Quaternion<Scalar>& q, const Vector<Scalar, 3>& acc_n)
{
    Vector<Scalar, 3> g_body = q.toRotationMatrix().transpose().col(2);
    return acc_n.cross(g_body);
}

// @cite mahony2008 -- Mahony et al., 2008, Eq. 48 (magnetic field error)
template <typename Scalar>
inline Vector<Scalar, 3> compute_magnetic_correction(const Eigen::Quaternion<Scalar>& q, const Vector<Scalar, 3>& mag_n)
{
    auto R = q.toRotationMatrix();
    Vector<Scalar, 3> h = R * mag_n;
    Vector<Scalar, 3> b{std::sqrt(h(0) * h(0) + h(1) * h(1)), Scalar{0}, h(2)};
    Vector<Scalar, 3> b_body = R.transpose() * b;
    return mag_n.cross(b_body);
}

}

#endif
