#ifndef HPP_GUARD_CTRLPP_ESTIMATION_MEKF_TYPES_H
#define HPP_GUARD_CTRLPP_ESTIMATION_MEKF_TYPES_H

#include "ctrlpp/types.h"

#include "ctrlpp/util/concepts.h"

#include <Eigen/Geometry>

#include <cmath>
#include <limits>
#include <cstddef>

namespace ctrlpp
{

template <typename M, typename Scalar, std::size_t NB, std::size_t NY>
concept mekf_measurement_model = requires(const M& m, const Eigen::Quaternion<Scalar>& q, const Vector<Scalar, NB>& b) {
    { m(q, b) } -> std::convertible_to<Vector<Scalar, NY>>;
};

template <typename M, typename Scalar, std::size_t NB, std::size_t NY>
concept differentiable_mekf_measurement = mekf_measurement_model<M, Scalar, NB, NY> && requires(const M& m, const Eigen::Quaternion<Scalar>& q, const Vector<Scalar, NB>& b) {
    { m.jacobian(q, b) } -> std::convertible_to<Matrix<Scalar, NY, 3 + NB>>;
};

enum class mekf_update_error
{
    non_finite_state,
    non_finite_covariance,
    non_finite_measurement,
    non_finite_result,
};

enum class mekf_health
{
    ok,
    non_finite_estimate,
};

template <ctrlpp_floating_scalar Scalar, std::size_t NB, std::size_t NY>
struct mekf_config
{
    static_assert(NB >= 3, "Bias dimension NB must be at least 3: the propagation subtracts the leading three bias elements from the gyro rate");
    static_assert(NY > 0, "Output dimension NY must be positive");
    static constexpr std::size_t NE = 3 + NB;
    Matrix<Scalar, NE, NE> Q{Matrix<Scalar, NE, NE>::Identity()};
    Matrix<Scalar, NY, NY> R{Matrix<Scalar, NY, NY>::Identity()};
    Eigen::Quaternion<Scalar> q0{Eigen::Quaternion<Scalar>::Identity()};
    Vector<Scalar, NB> b0{Vector<Scalar, NB>::Zero()};
    Matrix<Scalar, NE, NE> P0{Matrix<Scalar, NE, NE>::Identity()};
    Scalar dt{static_cast<Scalar>(0.01)};
    Scalar numerical_eps{std::cbrt(std::numeric_limits<Scalar>::epsilon())};
};

}

#endif
