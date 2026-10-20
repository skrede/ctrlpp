#ifndef HPP_GUARD_CTRLPP_ESTIMATION_SIGMA_POINTS_JULIER_SIGMA_POINTS_H
#define HPP_GUARD_CTRLPP_ESTIMATION_SIGMA_POINTS_JULIER_SIGMA_POINTS_H

/// @brief Symmetric sigma point generation (Julier-Uhlmann variant).
///
/// @cite julier2004 -- Julier & Uhlmann, "Unscented Filtering and Nonlinear Estimation", 2004

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/detail/covariance_ops.h"

#include "ctrlpp/estimation/estimation_types.h"
#include "ctrlpp/estimation/sigma_points/sigma_point_strategy.h"

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

public:
    static constexpr std::size_t num_points = 2 * NX + 1;
    using options_t = julier_options<Scalar>;

    /// @brief Construct from the default options, which lie inside the domain
    /// by construction and therefore need no validation.
    ///
    /// With the default zero kappa the sum n + kappa reduces to n = NX, which
    /// the class-scope assertion above requires to be positive.
    julier_sigma_points() : julier_sigma_points{unchecked_t{}, options_t{}} {}

    /// @brief Fallible factory. Validates the one exact domain condition the
    /// unscented transform imposes on kappa.
    ///
    /// The sum n + kappa is both the radicand of the offset scale
    /// gamma = sqrt(n + kappa) and the weight denominator, which fixes the
    /// admissible domain exactly: there is no tolerance here and no fitted
    /// constant. A NaN/Inf kappa, or a non-positive n + kappa, is rejected with
    /// filter_error::non_positive_scaling_radicand.
    ///
    /// A negative kappa with n + kappa > 0 lies inside the domain. Julier's own
    /// choice kappa = 3 - n is such a kappa for n > 3; it gives a negative center
    /// weight, which the method permits.
    ///
    /// @cite julier1997 -- Julier & Uhlmann, "A New Extension of the Kalman Filter to Nonlinear Systems", SPIE, 1997, Eq. 12
    static auto try_create(options_t opts = options_t{}) -> ctrlpp::expected<julier_sigma_points, filter_error>
    {
        if(!std::isfinite(opts.kappa) || !(n + opts.kappa > Scalar{0}))
            return ctrlpp::unexpected(filter_error::non_positive_scaling_radicand);
        return julier_sigma_points{unchecked_t{}, opts};
    }

    /// @brief Generate symmetric sigma point set with single tuning parameter kappa.
    ///
    /// @cite julier1997 -- Julier & Uhlmann, "A New Extension of the Kalman Filter to Nonlinear Systems", SPIE, 1997, Eq. 12
    /// @cite julier2004 -- Julier & Uhlmann, "Unscented Filtering and Nonlinear Estimation", 2004
    sigma_result<Scalar, NX, num_points> generate(const Vector<Scalar, NX>& x, const Matrix<Scalar, NX, NX>& P) const
    {
        sigma_result<Scalar, NX, num_points> result;

        Scalar gamma = std::sqrt(n + m_kappa);

        // Matrix square root factor S with S*S^T = P. The unpivoted Cholesky
        // factor is used so the reconstruction squares back to P exactly; a
        // non-positive-definite P is repaired to the nearest such matrix and
        // the event is surfaced through the result.
        auto sqrt_result = detail::covariance_sqrt<Scalar, nx>(P);
        const Eigen::Matrix<Scalar, nx, nx> S = sqrt_result.factor;
        result.spd_repaired = sqrt_result.repaired;

        // Center sigma point
        result.points[0] = x;

        // Offset sigma points
        for(std::size_t i = 0; i < NX; ++i)
        {
            Vector<Scalar, NX> offset = gamma * S.col(static_cast<int>(i));
            result.points[1 + i] = x + offset;
            result.points[1 + NX + i] = x - offset;
        }

        // Weights: Julier uses Wm = Wc
        Scalar denom = n + m_kappa;
        result.Wm[0] = m_kappa / denom;
        result.Wc[0] = m_kappa / denom;

        Scalar wi = Scalar{1} / (Scalar{2} * denom);
        for(std::size_t i = 1; i < num_points; ++i)
        {
            result.Wm[i] = wi;
            result.Wc[i] = wi;
        }

        return result;
    }

private:
    /// @brief Tag selecting the non-validating constructor reserved for
    /// `try_create` and for the in-domain default construction.
    struct unchecked_t
    {
        explicit unchecked_t() = default;
    };

    /// @brief Construct from options already known to be inside the domain.
    julier_sigma_points(unchecked_t, options_t opts) : m_kappa{opts.kappa} {}

    Scalar m_kappa;
};

}

#endif
