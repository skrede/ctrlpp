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

inline constexpr covariance_faults measurement_noise_faults{
        filter_error::non_finite_measurement_noise, filter_error::asymmetric_measurement_noise, filter_error::indefinite_measurement_noise};

inline constexpr covariance_faults initial_covariance_faults{
        filter_error::non_finite_initial_covariance, filter_error::asymmetric_initial_covariance, filter_error::indefinite_initial_covariance};

// M is written only when its asymmetry is drift, so an indefinite refusal leaves it symmetrized.
template <typename Scalar, int N>
ctrlpp::expected<psd_decomposition_result<Scalar, N>, filter_error> checked_covariance(Eigen::Matrix<Scalar, N, N>& M, const covariance_faults& faults)
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

// factor = 1 covers one working-precision normalization, at most (N / 2 + 2) u, plus this norm's own
// (N / 2 + 1) u: 7 u within N eps = 8 u for N = 4 (Higham, "Accuracy and Stability of Numerical
// Algorithms", 2nd ed., 2002, Sec. 3.1).
template <typename Scalar>
ctrlpp::expected<void, filter_error> checked_unit_quaternion(const Eigen::Quaternion<Scalar>& q, Scalar factor = Scalar{1})
{
    const Scalar norm = q.norm();
    if(!(norm > Scalar{0}) || !std::isfinite(norm))
        return ctrlpp::unexpected(filter_error::degenerate_quaternion);
    if(!(std::abs(norm - Scalar{1}) <= rounding_bound(4, Scalar{1}, factor)))
        return ctrlpp::unexpected(filter_error::non_unit_quaternion);
    return {};
}

// lambda_min > N eps lambda_max, written so that a NaN spectrum is singular.
template <typename Scalar, int N>
bool nonsingular_spectrum(const psd_decomposition_result<Scalar, N>& d)
{
    return d.eigenvalues.minCoeff() > static_cast<Scalar>(N) * std::numeric_limits<Scalar>::epsilon() * d.eigenvalues.maxCoeff();
}

// A refusal of R can follow a symmetrized Q; a caller that assigns nothing on refusal is unaffected.
template <typename Scalar, int NQ, int NR>
ctrlpp::expected<void, filter_error> validate_noise(Eigen::Matrix<Scalar, NQ, NQ>& Q, Eigen::Matrix<Scalar, NR, NR>& R)
{
    if(const auto process = checked_covariance(Q, process_noise_faults); !process)
        return ctrlpp::unexpected(process.error());
    if(const auto measurement = checked_covariance(R, measurement_noise_faults); !measurement)
        return ctrlpp::unexpected(measurement.error());
    return {};
}

template <typename Scalar, int N>
ctrlpp::expected<void, filter_error> validate_initial_covariance(Eigen::Matrix<Scalar, N, N>& P0)
{
    if(const auto prior = checked_covariance(P0, initial_covariance_faults); !prior)
        return ctrlpp::unexpected(prior.error());
    return {};
}

// The fields are independent, so the report order is their declaration order, not a claim about
// causality. On success Q, R and P0 hold their symmetrized forms.
template <typename Process, typename Measurement, typename State, typename Covariance>
ctrlpp::expected<void, filter_error> validate_filter_configuration(Process& Q, Measurement& R, const State& x0, Covariance& P0)
{
    if(const auto noise = validate_noise(Q, R); !noise)
        return ctrlpp::unexpected(noise.error());
    if(!x0.allFinite())
        return ctrlpp::unexpected(filter_error::non_finite_initial_state);
    return validate_initial_covariance(P0);
}

}

#endif
