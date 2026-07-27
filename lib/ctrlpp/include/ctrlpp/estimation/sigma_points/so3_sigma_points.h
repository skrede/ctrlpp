#ifndef HPP_GUARD_CTRLPP_ESTIMATION_SIGMA_POINTS_SO3_SIGMA_POINTS_H
#define HPP_GUARD_CTRLPP_ESTIMATION_SIGMA_POINTS_SO3_SIGMA_POINTS_H

/// @brief SO(3) manifold sigma point generation via tangent-space lifting.
///
/// @cite sola2018 -- Sola et al., "A micro Lie theory for state estimation in robotics", 2018
/// @cite hauberg2013 -- Hauberg et al., "Unscented Kalman Filtering on (Sub)Riemannian Manifolds", 2013

#include "ctrlpp/types.h"
#include "ctrlpp/config.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/lie/so3.h"

#include "ctrlpp/estimation/estimation_types.h"
#include "ctrlpp/estimation/sigma_points/merwe_sigma_points.h"
#include "ctrlpp/estimation/sigma_points/sigma_point_strategy.h"

#include <Eigen/Geometry>

#include <array>
#include <cstddef>
#include <concepts>
#include <utility>

namespace ctrlpp
{

template <typename Scalar, std::size_t NP>
struct manifold_sigma_result
{
    std::array<Eigen::Quaternion<Scalar>, NP> points;
    std::array<Scalar, NP> Wm;
    std::array<Scalar, NP> Wc;
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

public:
    static constexpr std::size_t num_points = merwe_sigma_points<Scalar, tangent_dim>::num_points;
    using options_t = merwe_options<Scalar>;

    /// @brief Construct from the default options, which the tangent-space
    /// strategy accepts by construction.
    so3_merwe_sigma_points() = default;

    /// @brief Fallible factory. Forwards the tangent-space strategy's parameter
    /// domain check unchanged, since the lifted points are generated from that
    /// strategy and inherit its weights.
    [[nodiscard]] static auto try_create(options_t opts = options_t{}) -> ctrlpp::expected<so3_merwe_sigma_points, filter_error>
    {
        auto inner = merwe_sigma_points<Scalar, tangent_dim>::try_create(opts);
        if(!inner)
            return ctrlpp::unexpected(inner.error());
        return so3_merwe_sigma_points{unchecked_t{}, std::move(*inner)};
    }

#if CTRLPP_HAS_EXCEPTIONS
    /// @brief Throwing convenience wrapper over `try_create`.
    ///
    /// Delegates to `try_create(opts).value()`. Compiled out when
    /// CTRLPP_HAS_EXCEPTIONS is 0; prefer `try_create` on exception-free builds.
    explicit so3_merwe_sigma_points(options_t opts) : so3_merwe_sigma_points{try_create(opts).value()} {}
#endif

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

        return result;
    }

private:
    /// @brief Tag selecting the non-validating constructor reserved for `try_create`.
    struct unchecked_t
    {
        explicit unchecked_t() = default;
    };

    /// @brief Construct from a tangent-space strategy already validated by `try_create`.
    so3_merwe_sigma_points(unchecked_t, merwe_sigma_points<Scalar, tangent_dim> inner) : m_inner{std::move(inner)} {}

    merwe_sigma_points<Scalar, tangent_dim> m_inner;
};

static_assert(manifold_sigma_point_strategy<so3_merwe_sigma_points<double>, double>);

}

#endif
