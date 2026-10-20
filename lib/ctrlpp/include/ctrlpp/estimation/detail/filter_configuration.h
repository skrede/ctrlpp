#ifndef HPP_GUARD_CTRLPP_ESTIMATION_DETAIL_FILTER_CONFIGURATION_H
#define HPP_GUARD_CTRLPP_ESTIMATION_DETAIL_FILTER_CONFIGURATION_H

#include "ctrlpp/estimation/estimation_types.h"

#include "ctrlpp/expected.h"

#include "ctrlpp/detail/covariance_ops.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <cmath>
#include <limits>

namespace ctrlpp::detail
{

struct covariance_faults
{
    filter_error non_finite;
    filter_error asymmetric;
    filter_error indefinite;
};

inline constexpr covariance_faults process_noise_faults{filter_error::non_finite_process_noise, filter_error::asymmetric_process_noise, filter_error::indefinite_process_noise};

inline constexpr covariance_faults measurement_noise_faults{filter_error::non_finite_measurement_noise, filter_error::asymmetric_measurement_noise, filter_error::indefinite_measurement_noise};

inline constexpr covariance_faults initial_covariance_faults{filter_error::non_finite_initial_covariance, filter_error::asymmetric_initial_covariance, filter_error::indefinite_initial_covariance};

/// @brief Refuse a candidate covariance that is non-finite, then one whose
/// asymmetry exceeds rounding drift, then one that is not positive
/// semidefinite; on success M holds its symmetrized form and the decomposition
/// of that matrix is handed back for factoring.
///
/// M is written only when its asymmetry is drift, so on an `indefinite`
/// refusal it already holds its symmetrized form.
template <typename Scalar, int N>
auto checked_covariance(Eigen::Matrix<Scalar, N, N>& M, const covariance_faults& faults) -> ctrlpp::expected<psd_decomposition_result<Scalar, N>, filter_error>
{
    if(!M.allFinite())
        return ctrlpp::unexpected(faults.non_finite);
    if(!symmetrize_within_drift(M))
        return ctrlpp::unexpected(faults.asymmetric);
    auto decomposition = psd_decomposition(M);
    if(!decomposition.psd)
        return ctrlpp::unexpected(faults.indefinite);
    return decomposition;
}

/// @brief Refuse an initial quaternion whose norm is zero or non-finite, then
/// one whose norm is not one within rounding.
///
/// The defect | ||q|| - 1 | is drift when it is within
/// `rounding_bound(N, 1, factor)` with N = 4 components. The count behind
/// `factor = 1`: a quaternion normalized once in working precision has a norm
/// defect of at most (N / 2 + 2) u (the N-term sum of squares, the square root
/// and the division), and this check's own evaluation of the norm adds at most
/// (N / 2 + 1) u (Higham, "Accuracy and Stability of Numerical Algorithms",
/// 2nd ed., 2002, Sec. 3.1); their sum, 7 u, is within N eps = 8 u. The caller
/// normalizes an admitted q, which removes the drift.
template <typename Scalar>
auto checked_unit_quaternion(const Eigen::Quaternion<Scalar>& q, Scalar factor = Scalar{1}) -> ctrlpp::expected<void, filter_error>
{
    const Scalar norm = q.norm();
    if(!(norm > Scalar{0}) || !std::isfinite(norm))
        return ctrlpp::unexpected(filter_error::degenerate_quaternion);
    if(!(std::abs(norm - Scalar{1}) <= rounding_bound(4, Scalar{1}, factor)))
        return ctrlpp::unexpected(filter_error::non_unit_quaternion);
    return {};
}

/// @brief Whether a decomposed positive semidefinite matrix is invertible:
/// lambda_min > N eps lambda_max, written so that a NaN spectrum is singular.
template <typename Scalar, int N>
auto nonsingular_spectrum(const psd_decomposition_result<Scalar, N>& d) -> bool
{
    return d.eigenvalues.minCoeff() > static_cast<Scalar>(N) * std::numeric_limits<Scalar>::epsilon() * d.eigenvalues.maxCoeff();
}

/// @brief Refuse a process noise Q, then a measurement noise R, that is
/// non-finite, structurally asymmetric or indefinite.
///
/// Q may be written by the first check and R by the second, so on success the
/// matrices the caller stores are the symmetrized ones. A refusal of R can
/// follow a projected Q; a caller that assigns nothing on refusal is unaffected.
template <typename Scalar, int NQ, int NR>
auto validate_noise(Eigen::Matrix<Scalar, NQ, NQ>& Q, Eigen::Matrix<Scalar, NR, NR>& R) -> ctrlpp::expected<void, filter_error>
{
    if(const auto process = checked_covariance(Q, process_noise_faults); !process)
        return ctrlpp::unexpected(process.error());
    if(const auto measurement = checked_covariance(R, measurement_noise_faults); !measurement)
        return ctrlpp::unexpected(measurement.error());
    return {};
}

/// @brief Refuse an initial covariance P0 that is non-finite, structurally
/// asymmetric or indefinite; on success P0 holds its symmetrized form.
template <typename Scalar, int N>
auto validate_initial_covariance(Eigen::Matrix<Scalar, N, N>& P0) -> ctrlpp::expected<void, filter_error>
{
    if(const auto prior = checked_covariance(P0, initial_covariance_faults); !prior)
        return ctrlpp::unexpected(prior.error());
    return {};
}

/// @brief Validate the noise and initial-condition fields that `kalman_config`,
/// `ekf_config` and `ukf_config` declare identically.
///
/// The fields are independent of one another, unlike a step's carried estimate
/// and its measurement, so the report order is the order the aggregates
/// declare, not a claim about causality. Each matrix is refused when it is
/// non-finite, then when it is structurally asymmetric, then when it is
/// indefinite; a singular positive semidefinite matrix is accepted. On success
/// Q, R and P0 hold their symmetrized forms, which are what the caller stores.
template <typename Process, typename Measurement, typename State, typename Covariance>
auto validate_filter_configuration(Process& Q, Measurement& R, const State& x0, Covariance& P0) -> ctrlpp::expected<void, filter_error>
{
    if(const auto noise = validate_noise(Q, R); !noise)
        return ctrlpp::unexpected(noise.error());
    if(!x0.allFinite())
        return ctrlpp::unexpected(filter_error::non_finite_initial_state);
    return validate_initial_covariance(P0);
}

}

#endif
