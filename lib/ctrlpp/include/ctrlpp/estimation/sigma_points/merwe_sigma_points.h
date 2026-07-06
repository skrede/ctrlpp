#ifndef HPP_GUARD_CTRLPP_ESTIMATION_SIGMA_POINTS_MERWE_SIGMA_POINTS_H
#define HPP_GUARD_CTRLPP_ESTIMATION_SIGMA_POINTS_MERWE_SIGMA_POINTS_H

/// @brief Scaled sigma point generation (Van der Merwe variant).
///
/// @cite vandermerwe2004 -- Van der Merwe, "Sigma-Point Kalman Filters", PhD thesis, 2004

#include "ctrlpp/types.h"

#include "ctrlpp/detail/covariance_ops.h"

#include "ctrlpp/estimation/sigma_points/sigma_point_strategy.h"

#include <cmath>
#include <cstddef>

namespace ctrlpp
{

template <typename Scalar>
struct merwe_options
{
    Scalar alpha{Scalar{1e-3}};
    Scalar beta{Scalar{2}};
    Scalar kappa{Scalar{0}};
};

template <typename Scalar, std::size_t NX>
class merwe_sigma_points
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr Scalar n = static_cast<Scalar>(NX);

public:
    static constexpr std::size_t num_points = 2 * NX + 1;
    using options_t = merwe_options<Scalar>;

    explicit merwe_sigma_points(options_t opts = options_t{}) : m_alpha{opts.alpha}, m_beta{opts.beta}, m_kappa{opts.kappa} {}

    /// @brief Generate scaled symmetric sigma point set with mean/covariance weights.
    ///
    /// @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001, Eq. 15
    /// @cite vandermerwe2004 -- Van der Merwe, "Sigma-Point Kalman Filters", PhD thesis, 2004
    sigma_result<Scalar, NX, num_points> generate(const Vector<Scalar, NX>& x, const Matrix<Scalar, NX, NX>& P) const
    {
        sigma_result<Scalar, NX, num_points> result;

        Scalar lambda = m_alpha * m_alpha * (n + m_kappa) - n;
        Scalar gamma = std::sqrt(n + lambda);

        // Matrix square root factor S with S*S^T = P. The unpivoted Cholesky
        // factor is used so the reconstruction squares back to P exactly; a
        // non-positive-definite P is repaired to the nearest such matrix and
        // the event is surfaced through the result.
        auto sqrt_result = detail::covariance_sqrt<Scalar, nx>(P);
        const Eigen::Matrix<Scalar, nx, nx> S = sqrt_result.factor;
        result.spd_repaired = sqrt_result.repaired;

        // Center sigma point
        result.points[0] = x;

        // Offset sigma points: x +/- gamma * columns of S
        for(std::size_t i = 0; i < NX; ++i)
        {
            Vector<Scalar, NX> offset = gamma * S.col(static_cast<int>(i));
            result.points[1 + i] = x + offset;
            result.points[1 + NX + i] = x - offset;
        }

        // Weights
        Scalar denom = n + lambda;
        result.Wm[0] = lambda / denom;
        result.Wc[0] = lambda / denom + (Scalar{1} - m_alpha * m_alpha + m_beta);

        Scalar wi = Scalar{1} / (Scalar{2} * denom);
        for(std::size_t i = 1; i < num_points; ++i)
        {
            result.Wm[i] = wi;
            result.Wc[i] = wi;
        }

        return result;
    }

private:
    Scalar m_alpha;
    Scalar m_beta;
    Scalar m_kappa;
};

}

#endif
