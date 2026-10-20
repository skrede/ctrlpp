#ifndef HPP_GUARD_CTRLPP_ESTIMATION_ESTIMATION_TYPES_H
#define HPP_GUARD_CTRLPP_ESTIMATION_ESTIMATION_TYPES_H

/// @brief Public types shared by the estimation filters.
///
/// `filter_error` enumerates the structured failure modes of the fallible
/// factories in the estimation module: the quaternion attitude filters (`mekf`,
/// `manifold_ukf`, `complementary_filter`), the sigma-point strategies, the
/// filters that build one from an options aggregate, and the configuration
/// validation of every filter that takes a Q, R or P0, which refuses a matrix
/// that is non-finite, structurally asymmetric or not positive semidefinite.
/// It forms their `ctrlpp::expected<T, filter_error>` construction contract.

namespace ctrlpp
{

/// @brief Structured failure modes for the estimation-module factories.
///
///  * degenerate_quaternion         : the initial quaternion q0 has zero or
///                                    non-finite norm, so normalizing it would
///                                    produce NaN and silently poison the filter
///                                    state at construction. It is reported
///                                    ahead of the unit-norm check, which a
///                                    finite nonzero q0 is judged by next.
///  * non_positive_sigma_spread     : the sigma-point spread parameter is not
///                                    finite and strictly positive. It is the
///                                    divisor of the sigma-point weight
///                                    denominator, so this is an exact domain
///                                    condition and not a tuning preference: a
///                                    zero spread makes the weights infinite,
///                                    and a negative one makes them finite but
///                                    wrong, since only its square is used.
///  * non_positive_scaling_radicand : the sum of the state dimension and the
///                                    secondary scaling parameter is not finite
///                                    and strictly positive. It sits under the
///                                    square root that scales the sigma-point
///                                    offsets and inside the same weight
///                                    denominator, so a non-positive sum yields
///                                    non-finite offsets or non-finite weights.
///  * non_finite_process_noise      : the process-noise matrix Q has a
///                                    non-finite entry. Q is added to the
///                                    propagated covariance every prediction, so
///                                    an infinite entry makes the predicted
///                                    covariance infinite, the gain a ratio of
///                                    infinities and therefore non-finite, and
///                                    the corrected estimate non-finite with it.
///  * non_finite_measurement_noise  : the measurement-noise matrix R has a
///                                    non-finite entry. R is added to the
///                                    innovation covariance, which is the matrix
///                                    the gain solve is posed against, so a
///                                    non-finite entry makes that solve
///                                    meaningless for every measurement.
///  * non_finite_initial_state      : the initial state vector x0 has a
///                                    non-finite component. It seeds the carried
///                                    estimate, which is the filter's memory, so
///                                    every later estimate is formed from it.
///  * non_finite_initial_covariance : the initial covariance matrix P0 has a
///                                    non-finite entry. It seeds the covariance
///                                    recursion, and the recursion has no
///                                    mechanism that could return a non-finite
///                                    covariance to a finite one.
///  * indefinite_process_noise      : Q has an eigenvalue below -floor, with
///                                    floor = max(N eps max|Q_ij|, N denorm_min):
///                                    the symmetric eigensolver's backward-error
///                                    order, the convention
///                                    `detail::covariance_sqrt` cites, kept from
///                                    underflowing to zero at subnormal scale.
///                                    Such a matrix is not a covariance: there
///                                    is no Gaussian with it to propagate or to
///                                    sample the process noise from.
///  * indefinite_measurement_noise  : R has an eigenvalue below the same floor,
///                                    so there is no Gaussian density with it to
///                                    score a measurement against.
///  * singular_measurement_noise    : lambda_min(R) <= NY eps lambda_max(R), the
///                                    reciprocal-pivot convention the continuous
///                                    LQR gain applies to its input weight. R
///                                    then has no inverse and no finite
///                                    log-determinant, both of which a filter
///                                    that evaluates a Gaussian density needs.
///  * indefinite_initial_covariance : P0 has an eigenvalue below the same floor,
///                                    so there is no Gaussian prior with it.
///  * asymmetric_process_noise      : max|Q_ij - Q_ji| exceeds the rounding
///                                    bound max(N eps max|Q_ij|, N denorm_min).
///                                    A covariance is symmetric by definition,
///                                    and an asymmetry beyond the rounding of
///                                    the products that formed Q is not drift
///                                    but a different matrix in each triangle.
///  * asymmetric_measurement_noise  : R's asymmetry exceeds the same bound.
///  * asymmetric_initial_covariance : P0's asymmetry exceeds the same bound.
///  * non_unit_quaternion           : | ||q0|| - 1 | exceeds the rounding bound
///                                    max(4 eps, 4 denorm_min). A rotation
///                                    quaternion is unit by definition, and a
///                                    norm defect beyond the rounding of a
///                                    normalization is not drift but a
///                                    quaternion that names no attitude. A
///                                    defect within the bound is normalized
///                                    away at construction.
///  * non_finite_initial_bias       : the initial bias b0 of the multiplicative
///                                    EKF has a non-finite component. It seeds
///                                    the bias estimate, which is part of the
///                                    carried state, so every later attitude
///                                    is formed from it.
///
/// Every filter that takes a Q, R or P0 refuses one that is non-finite,
/// structurally asymmetric or indefinite, in that order for each matrix: a
/// matrix that is not symmetric positive semidefinite is not a covariance,
/// whatever the filter does with it. An asymmetry within the rounding bound is
/// drift: it is accepted, the matrix the filter stores is its symmetrized form
/// (each differing pair replaced by its midpoint), and the definiteness verdict
/// runs on that symmetrized matrix. A singular positive semidefinite matrix is
/// accepted: a point-mass prior, a rank-one process noise and a noiseless
/// measurement channel are legitimate configurations, and the Kalman-family
/// filters only add and multiply Q, R and P0, inverting the innovation
/// covariance and never R alone. An ill-conditioned configuration, such as a
/// covariance with entries many orders of magnitude apart, is accepted as well:
/// conditioning is a numerical-behavior question, not a domain condition.
///
/// The particle filter additionally refuses a singular R, because it scores
/// every particle with R's inverse and log-determinant.
enum class filter_error
{
    degenerate_quaternion,
    non_positive_sigma_spread,
    non_positive_scaling_radicand,
    non_finite_process_noise,
    non_finite_measurement_noise,
    non_finite_initial_state,
    non_finite_initial_covariance,
    indefinite_process_noise,
    indefinite_measurement_noise,
    singular_measurement_noise,
    indefinite_initial_covariance,
    asymmetric_process_noise,
    asymmetric_measurement_noise,
    asymmetric_initial_covariance,
    non_unit_quaternion,
    non_finite_initial_bias,
};

}

#endif
