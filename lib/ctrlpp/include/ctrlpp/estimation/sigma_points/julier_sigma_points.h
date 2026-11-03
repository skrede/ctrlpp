#ifndef HPP_GUARD_CTRLPP_ESTIMATION_SIGMA_POINTS_JULIER_SIGMA_POINTS_H
#define HPP_GUARD_CTRLPP_ESTIMATION_SIGMA_POINTS_JULIER_SIGMA_POINTS_H

// @cite julier2004 -- Julier & Uhlmann, "Unscented Filtering and Nonlinear Estimation", 2004

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
struct julier_options
{
    Scalar kappa{Scalar{0}};
};

template <typename Scalar, std::size_t NX>
class julier_sigma_points
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
    using options_t = julier_options<Scalar>;

    // The default zero kappa reduces n + kappa to NX, which the class assertion requires to be positive,
    // so the default options need no validation.
    julier_sigma_points()
        : julier_sigma_points{unchecked_t{}, options_t{}}
    {
    }

    // n + kappa is both the radicand of gamma = sqrt(n + kappa) and the weight denominator, so its
    // positivity is the exact domain, with no tolerance.
    // @cite julier1997 -- Julier & Uhlmann, "A New Extension of the Kalman Filter to Nonlinear Systems", SPIE, 1997, Eq. 12
    static ctrlpp::expected<julier_sigma_points, filter_error> try_create(options_t opts = options_t{})
    {
        if(!std::isfinite(opts.kappa) || !(n + opts.kappa > Scalar{0}))
            return ctrlpp::unexpected(filter_error::non_positive_scaling_radicand);
        return julier_sigma_points{unchecked_t{}, opts};
    }

    // The unpivoted Cholesky factor squares back to P exactly; a P that is not positive definite is
    // repaired to the nearest such matrix and the repair is reported through the result.
    // @cite julier1997 -- Julier & Uhlmann, "A New Extension of the Kalman Filter to Nonlinear Systems", SPIE, 1997, Eq. 12
    // @cite julier2004 -- Julier & Uhlmann, "Unscented Filtering and Nonlinear Estimation", 2004
    sigma_result<Scalar, NX, num_points> generate(const Vector<Scalar, NX>& x, const Matrix<Scalar, NX, NX>& P) const
    {
        sigma_result<Scalar, NX, num_points> result;
        Scalar gamma = std::sqrt(n + m_kappa);

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

        Scalar denom = n + m_kappa;
        result.Wm[0] = m_kappa / denom;
        result.Wc[0] = m_kappa / denom;
        Scalar wi = Scalar{1} / (Scalar{2} * denom);
        fill_offset_weights(result, wi);
        return result;
    }

private:
    Scalar m_kappa;

    julier_sigma_points(unchecked_t, options_t opts)
        : m_kappa{opts.kappa}
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
