#ifndef HPP_GUARD_CTRLPP_ESTIMATION_MANIFOLD_UKF_TYPES_H
#define HPP_GUARD_CTRLPP_ESTIMATION_MANIFOLD_UKF_TYPES_H

#include "ctrlpp/types.h"

#include "ctrlpp/util/concepts.h"

#include <Eigen/Geometry>

#include <cstddef>
#include <cstdint>
#include <concepts>

namespace ctrlpp
{

template <typename M, typename Scalar, std::size_t NY>
concept manifold_ukf_measurement_model = requires(const M& m, const Eigen::Quaternion<Scalar>& q) {
    { m(q) } -> std::convertible_to<Vector<Scalar, NY>>;
};

template <typename D, typename Scalar>
concept manifold_ukf_dynamics_model = requires(const D& d, const Eigen::Quaternion<Scalar>& q, const Vector<Scalar, 3>& u) {
    { d(q, u) } -> std::convertible_to<Eigen::Quaternion<Scalar>>;
};

enum class manifold_ukf_update_error
{
    non_finite_state,
    non_finite_covariance,
    non_finite_measurement,
    non_finite_result,
};

// Ordered by severity, which note_prediction relies on; the status latches and never returns to a lower rung.
enum class manifold_ukf_health
{
    ok,
    covariance_repaired,
    mean_not_converged,
    non_finite_estimate
};

template <ctrlpp_floating_scalar Scalar, std::size_t NY>
struct manifold_ukf_config
{
    static_assert(NY > 0, "Output dimension NY must be positive");
    Matrix<Scalar, 3, 3> Q{Matrix<Scalar, 3, 3>::Identity()};
    Matrix<Scalar, NY, NY> R{Matrix<Scalar, NY, NY>::Identity()};
    Eigen::Quaternion<Scalar> q0{Eigen::Quaternion<Scalar>::Identity()};
    Matrix<Scalar, 3, 3> P0{Matrix<Scalar, 3, 3>::Identity()};
    std::int32_t geodesic_mean_max_iter{30};
    Scalar geodesic_mean_tol{static_cast<Scalar>(1e-9)};
};

}

#endif
