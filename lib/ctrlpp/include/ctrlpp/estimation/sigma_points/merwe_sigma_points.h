#ifndef HPP_GUARD_CTRLPP_ESTIMATION_SIGMA_POINTS_MERWE_SIGMA_POINTS_H
#define HPP_GUARD_CTRLPP_ESTIMATION_SIGMA_POINTS_MERWE_SIGMA_POINTS_H

/// @brief Scaled sigma point generation (Van der Merwe variant).
///
/// @cite vandermerwe2004 -- Van der Merwe, "Sigma-Point Kalman Filters", PhD thesis, 2004

#include "ctrlpp/types.h"
#include "ctrlpp/config.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/detail/covariance_ops.h"

#include "ctrlpp/estimation/estimation_types.h"
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
    static_assert(NX > 0, "State dimension NX must be positive");

    static constexpr int nx = static_cast<int>(NX);
    static constexpr Scalar n = static_cast<Scalar>(NX);

public:
    static constexpr std::size_t num_points = 2 * NX + 1;
    using options_t = merwe_options<Scalar>;

    /// @brief Construct from the default options, which lie inside the domain
    /// by construction and therefore need no validation.
    ///
    /// The default spread is finite and strictly positive, and with the default
    /// zero kappa the dimension sum n + kappa reduces to n = NX, which the
    /// class-scope assertion above requires to be positive.
    merwe_sigma_points() : merwe_sigma_points{unchecked_t{}, options_t{}} {}

    /// @brief Fallible factory. Validates the two exact domain conditions the
    /// scaled unscented transform imposes on its parameters.
    ///
    /// The scaling term is lambda = alpha^2 (n + kappa) - n, so the weight
    /// denominator n + lambda and the squared sigma-point offset scale
    /// gamma^2 = n + lambda both collapse to alpha^2 (n + kappa). The divisor
    /// and the radicand are therefore the same expression, which fixes the
    /// admissible domain exactly: there is no tolerance here and no fitted
    /// constant. Rejections, checked in order:
    ///  * NaN/Inf or non-positive alpha -> filter_error::non_positive_sigma_spread
    ///  * NaN/Inf kappa, or non-positive n + kappa -> filter_error::non_positive_scaling_radicand
    ///
    /// A negative alpha is rejected alongside a zero one: only alpha^2 reaches
    /// the denominator, so a negative spread produces finite but wrong weights,
    /// which no finiteness check on the generated set would catch.
    ///
    /// Nothing else is validated. beta enters only the additive prior-kurtosis
    /// term of the first covariance weight and carries no domain restriction of
    /// this kind. Choosing good default values across dimension, scale, and
    /// scalar tier is a separate question from admissibility and is not decided
    /// here.
    ///
    /// @cite wan2001 -- Wan & van der Merwe, "The Unscented Kalman Filter", 2001, Eq. 15
    [[nodiscard]] static auto try_create(options_t opts = options_t{}) -> ctrlpp::expected<merwe_sigma_points, filter_error>
    {
        if(!std::isfinite(opts.alpha) || opts.alpha <= Scalar{0})
            return ctrlpp::unexpected(filter_error::non_positive_sigma_spread);
        if(!std::isfinite(opts.kappa) || !(n + opts.kappa > Scalar{0}))
            return ctrlpp::unexpected(filter_error::non_positive_scaling_radicand);
        return merwe_sigma_points{unchecked_t{}, opts};
    }

#if CTRLPP_HAS_EXCEPTIONS
    /// @brief Throwing convenience wrapper over `try_create`.
    ///
    /// Delegates to `try_create(opts).value()`, so an out-of-domain parameter
    /// set throws the value() exception of `ctrlpp::expected`. Compiled out
    /// when CTRLPP_HAS_EXCEPTIONS is 0; prefer `try_create` on exception-free
    /// builds, where it is the only construction path that takes options.
    explicit merwe_sigma_points(options_t opts) : merwe_sigma_points{try_create(opts).value()} {}
#endif

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
    /// @brief Tag selecting the non-validating constructor reserved for
    /// `try_create` and for the in-domain default construction.
    struct unchecked_t
    {
        explicit unchecked_t() = default;
    };

    /// @brief Construct from options already known to be inside the domain.
    merwe_sigma_points(unchecked_t, options_t opts) : m_alpha{opts.alpha}, m_beta{opts.beta}, m_kappa{opts.kappa} {}

    Scalar m_alpha;
    Scalar m_beta;
    Scalar m_kappa;
};

}

#endif
