#ifndef HPP_GUARD_CTRLPP_ESTIMATION_DETAIL_GEODESIC_MEAN_H
#define HPP_GUARD_CTRLPP_ESTIMATION_DETAIL_GEODESIC_MEAN_H

#include "ctrlpp/types.h"

#include "ctrlpp/lie/so3.h"

#include <Eigen/Geometry>

#include <array>
#include <cstddef>
#include <cstdint>

namespace ctrlpp::detail
{

template <typename Scalar>
struct geodesic_mean_result
{
    Eigen::Quaternion<Scalar> mean;
    bool converged;
    Scalar residual;
};

// @cite hauberg2013 -- Hauberg et al., "Unscented Kalman Filtering on (Sub)Riemannian Manifolds", 2013, Alg. 1
template <typename Scalar, std::size_t NP>
geodesic_mean_result<Scalar> geodesic_mean_impl(const std::array<Eigen::Quaternion<Scalar>, NP>& qs, const std::array<Scalar, NP>& Wm, std::int32_t max_iter, Scalar tol)
{
    auto q_mean = qs[0];
    Scalar residual = Scalar{0};
    bool converged = false;
    for(std::int32_t iter = 0; iter < max_iter; ++iter)
    {
        Vector<Scalar, 3> eps = Vector<Scalar, 3>::Zero();
        for(std::size_t i = 0; i < NP; ++i)
        {
            Eigen::Quaternion<Scalar> qi = qs[i];
            if(qi.dot(q_mean) < Scalar{0})
            {
                qi.w() = -qi.w();
                qi.x() = -qi.x();
                qi.y() = -qi.y();
                qi.z() = -qi.z();
            }
            eps += Wm[i] * so3::log(q_mean.conjugate() * qi);
        }
        residual = eps.norm();
        if(residual < tol)
        {
            converged = true;
            break;
        }
        q_mean = (q_mean * so3::exp(eps)).normalized();
    }
    return {q_mean, converged, residual};
}

}

#endif
