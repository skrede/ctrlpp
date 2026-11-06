#ifndef HPP_GUARD_CTRLPP_ESTIMATION_SIGMA_POINTS_MERWE_SIGMA_POINTS_H
#define HPP_GUARD_CTRLPP_ESTIMATION_SIGMA_POINTS_MERWE_SIGMA_POINTS_H

// @cite vandermerwe2004 -- Van der Merwe, "Sigma-Point Kalman Filters", PhD thesis, 2004

#include "ctrlpp/estimation/sigma_points/sigma_point_strategy.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/detail/covariance_ops.h"

#include "ctrlpp/estimation/estimation_types.h"

#include <cmath>
#include <cstddef>

namespace ctrlpp
{

template <typename Scalar>
struct merwe_options
{
    Scalar alpha{static_cast<Scalar>(1e-3)};
    Scalar beta{Scalar{2}};
    Scalar kappa{Scalar{0}};
};

template <typename Scalar, std::size_t NX>
class merwe_sigma_points
{
    static_assert(NX > 0, "State dimension NX must be positive");

    static constexpr int nx = static_cast<int>(NX);
    static constexpr Scalar n = static_cast<Scalar>(NX);

    struct unchecked_t
    {
        explicit unchecked_t() = default;
    };

public:
    static constexpr std::size_t num_points = 2 * NX + 1;
    using options_t = merwe_options<Scalar>;

    // The default spread is positive and the default zero kappa reduces n + kappa to NX, which the class
    // assertion requires to be positive, so the default options need no validation.
    merwe_sigma_points()
        : merwe_sigma_points{unchecked_t{}, options_t{}}
    {
    }

    // With lambda = alpha^2 (n + kappa) - n, the weight denominator and the squared offset scale both equal
    // alpha^2 (n + kappa), so its positivity is the exact domain. A negative alpha is refused as well: only
    // alpha^2 reaches the weights, so it would give finite but wrong ones.
    // @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001, Eq. 7.30
    static ctrlpp::expected<merwe_sigma_points, filter_error> try_create(options_t opts = options_t{})
    {
        if(!std::isfinite(opts.alpha) || opts.alpha <= Scalar{0})
            return ctrlpp::unexpected(filter_error::non_positive_sigma_spread);
        if(!std::isfinite(opts.kappa) || !(n + opts.kappa > Scalar{0}))
            return ctrlpp::unexpected(filter_error::non_positive_scaling_radicand);
        return merwe_sigma_points{unchecked_t{}, opts};
    }

    // The unpivoted Cholesky factor squares back to P exactly; a P that is not positive definite is
    // repaired to the nearest such matrix and the repair is reported through the result.
    // @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001, Eq. 7.30, 7.34
    // @cite vandermerwe2004 -- Van der Merwe, "Sigma-Point Kalman Filters", PhD thesis, 2004
    sigma_result<Scalar, NX, num_points> generate(const Vector<Scalar, NX>& x, const Matrix<Scalar, NX, NX>& P) const
    {
        sigma_result<Scalar, NX, num_points> result;
        Scalar lambda = m_alpha * m_alpha * (n + m_kappa) - n;
        Scalar gamma = std::sqrt(n + lambda);

        auto sqrt_result = detail::covariance_sqrt<Scalar, nx>(P);
        const Eigen::Matrix<Scalar, nx, nx> S = sqrt_result.factor;
        result.spd_repaired = sqrt_result.repaired;
        result.factor = sqrt_result.factor;

        result.points[0] = x;
        for(std::size_t i = 0; i < NX; ++i)
        {
            Vector<Scalar, NX> offset = gamma * S.col(static_cast<Eigen::Index>(i));
            result.points[1 + i] = x + offset;
            result.points[1 + NX + i] = x - offset;
        }

        Scalar denom = n + lambda;
        result.Wm[0] = lambda / denom;
        result.Wc[0] = lambda / denom + (Scalar{1} - m_alpha * m_alpha + m_beta);
        Scalar wi = Scalar{1} / (Scalar{2} * denom);
        fill_offset_weights(result, wi);
        return result;
    }

private:
    Scalar m_alpha;
    Scalar m_beta;
    Scalar m_kappa;

    merwe_sigma_points(unchecked_t, options_t opts)
        : m_alpha{opts.alpha}
        , m_beta{opts.beta}
        , m_kappa{opts.kappa}
    {
    }

    static void fill_offset_weights(sigma_result<Scalar, NX, num_points>& result, Scalar wi)
    {
        for(std::size_t i = 1; i < num_points; ++i)
        {
            result.Wm[i] = wi;
            result.Wc[i] = wi;
        }
    }
};

}

#endif
