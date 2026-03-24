#ifndef HPP_GUARD_CTRLPP_ESTIMATION_SIGMA_POINTS_SO3_SIGMA_POINTS_H
#define HPP_GUARD_CTRLPP_ESTIMATION_SIGMA_POINTS_SO3_SIGMA_POINTS_H

#include "ctrlpp/lie/so3.h"
#include "ctrlpp/types.h"

#include "ctrlpp/estimation/sigma_points/merwe_sigma_points.h"
#include "ctrlpp/estimation/sigma_points/sigma_point_strategy.h"

#include <Eigen/Geometry>

#include <array>
#include <cstddef>
#include <concepts>

namespace ctrlpp {

template <typename Scalar, std::size_t NP>
struct manifold_sigma_result
{
    std::array<Eigen::Quaternion<Scalar>, NP> points;
    std::array<Scalar, NP> Wm;
    std::array<Scalar, NP> Wc;
};

template <typename S, typename Scalar>
concept manifold_sigma_point_strategy =
    requires
    {
        { S::num_points } -> std::convertible_to<std::size_t>;
        typename S::options_t;
    } &&
    requires(const S &s, const Eigen::Quaternion<Scalar> &q, const Matrix<Scalar, 3, 3> &P)
    {
        { s.generate(q, P) } -> std::convertible_to<manifold_sigma_result<Scalar, S::num_points>>;
    };

template <typename Scalar>
class so3_merwe_sigma_points
{
    static constexpr std::size_t tangent_dim = 3;

public:
    static constexpr std::size_t num_points = merwe_sigma_points<Scalar, tangent_dim>::num_points;
    using options_t = merwe_options<Scalar>;

    explicit so3_merwe_sigma_points(options_t opts = options_t{})
        : m_inner{opts}
    {
    }

    manifold_sigma_result<Scalar, num_points> generate(const Eigen::Quaternion<Scalar> &q_mean, const Matrix<Scalar, 3, 3> &P) const

    {
        auto tangent = m_inner.generate(Vector<Scalar, tangent_dim>::Zero(), P);

        manifold_sigma_result<Scalar, num_points> result;

        for(std::size_t i = 0; i < num_points; ++i)
        {
            result.points[i] = so3::compose(q_mean, so3::exp(tangent.points[i]));
            result.points[i].normalize();
        }

        result.Wm = tangent.Wm;
        result.Wc = tangent.Wc;

        return result;
    }

private:
    merwe_sigma_points<Scalar, tangent_dim> m_inner;
};

static_assert(manifold_sigma_point_strategy<so3_merwe_sigma_points<double>, double>);

}

#endif
