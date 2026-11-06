#ifndef HPP_GUARD_CTRLPP_ESTIMATION_SIGMA_POINTS_SO3_SIGMA_POINTS_H
#define HPP_GUARD_CTRLPP_ESTIMATION_SIGMA_POINTS_SO3_SIGMA_POINTS_H

// @cite sola2018 -- Sola et al., "A micro Lie theory for state estimation in robotics", 2018
// @cite hauberg2013 -- Hauberg, Lauze & Pedersen, "Unscented Kalman Filtering on Riemannian Manifolds", 2013

#include "ctrlpp/estimation/sigma_points/merwe_sigma_points.h"
#include "ctrlpp/estimation/sigma_points/sigma_point_strategy.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/lie/so3.h"

#include "ctrlpp/estimation/estimation_types.h"

#include <Eigen/Geometry>

#include <array>
#include <cstddef>
#include <utility>
#include <concepts>

namespace ctrlpp
{

template <typename Scalar, std::size_t NP>
struct manifold_sigma_result
{
    manifold_sigma_result()
        : spd_repaired{false}
    {
    }

    std::array<Eigen::Quaternion<Scalar>, NP> points;
    std::array<Scalar, NP> Wm;
    std::array<Scalar, NP> Wc;
    bool spd_repaired;
    // S with S*S^T equal to the tangent covariance the points were built from, the repaired one after a
    // repair. Read only when spd_repaired is set, so a strategy that sets the flag must fill it.
    Matrix<Scalar, 3, 3> factor;
};

template <typename S, typename Scalar>
concept manifold_sigma_point_strategy = requires {
    { S::num_points } -> std::convertible_to<std::size_t>;
    typename S::options_t;
} && requires(const S& s, const Eigen::Quaternion<Scalar>& q, const Matrix<Scalar, 3, 3>& P) {
    { s.generate(q, P) } -> std::convertible_to<manifold_sigma_result<Scalar, S::num_points>>;
};

template <typename Scalar>
class so3_merwe_sigma_points
{
    static constexpr std::size_t tangent_dim = 3;

    struct unchecked_t
    {
        explicit unchecked_t() = default;
    };

public:
    static constexpr std::size_t num_points = merwe_sigma_points<Scalar, tangent_dim>::num_points;
    using options_t = merwe_options<Scalar>;

    so3_merwe_sigma_points() = default;

    // The lifted points inherit the tangent-space strategy's weights, so its domain check is forwarded unchanged.
    static ctrlpp::expected<so3_merwe_sigma_points, filter_error> try_create(options_t opts = options_t{})
    {
        auto inner = merwe_sigma_points<Scalar, tangent_dim>::try_create(opts);
        if(!inner)
            return ctrlpp::unexpected(inner.error());
        return so3_merwe_sigma_points{unchecked_t{}, std::move(*inner)};
    }

    manifold_sigma_result<Scalar, num_points> generate(const Eigen::Quaternion<Scalar>& q_mean, const Matrix<Scalar, 3, 3>& P) const
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
        result.spd_repaired = tangent.spd_repaired;
        result.factor = tangent.factor;

        return result;
    }

private:
    merwe_sigma_points<Scalar, tangent_dim> m_inner;

    so3_merwe_sigma_points(unchecked_t, merwe_sigma_points<Scalar, tangent_dim> inner)
        : m_inner{std::move(inner)}
    {
    }
};

static_assert(manifold_sigma_point_strategy<so3_merwe_sigma_points<double>, double>);

}

#endif
