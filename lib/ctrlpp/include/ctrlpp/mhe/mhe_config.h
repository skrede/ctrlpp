#ifndef HPP_GUARD_CTRLPP_MHE_MHE_CONFIG_H
#define HPP_GUARD_CTRLPP_MHE_MHE_CONFIG_H

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/estimation/estimation_types.h"

#include "ctrlpp/mpc/nmpc_config.h"

#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <variant>
#include <optional>
#include <functional>

namespace ctrlpp
{

enum class moving_horizon_configuration_error : std::uint8_t
{
    non_invertible_process_noise,
    non_invertible_measurement_noise,
    non_invertible_initial_covariance,
    non_positive_arrival_cost_weight,
    invalid_state_bounds,
    invalid_residual_bound,
    non_positive_soft_penalty,
    non_positive_numerical_eps,
    invalid_path_penalty,
    non_finite_path_constraint,
};

using moving_horizon_construction_error = std::variant<filter_error, moving_horizon_configuration_error>;

namespace detail
{

template <typename MatrixType>
ctrlpp::expected<MatrixType, moving_horizon_configuration_error> finite_full_piv_inverse(const MatrixType& matrix, moving_horizon_configuration_error error)
{
    auto factor = matrix.fullPivLu();
    if(!factor.isInvertible())
        return ctrlpp::unexpected(error);

    MatrixType inverse = factor.solve(MatrixType::Identity());
    if(!inverse.allFinite())
        return ctrlpp::unexpected(error);
    return inverse;
}

template <typename Config>
ctrlpp::expected<void, moving_horizon_configuration_error> validate_moving_horizon_options(const Config& config)
{
    using scalar_type = typename Config::scalar_type;

    if(!std::isfinite(config.arrival_cost_weight) || !(config.arrival_cost_weight > scalar_type{0}))
        return ctrlpp::unexpected(moving_horizon_configuration_error::non_positive_arrival_cost_weight);
    if((config.x_min && !config.x_min->allFinite()) || (config.x_max && !config.x_max->allFinite()) ||
            (config.x_min && config.x_max && !((*config.x_min).array() <= (*config.x_max).array()).all()))
        return ctrlpp::unexpected(moving_horizon_configuration_error::invalid_state_bounds);
    if(config.residual_bound && (!config.residual_bound->allFinite() || !(config.residual_bound->array() >= 0).all()))
        return ctrlpp::unexpected(moving_horizon_configuration_error::invalid_residual_bound);
    if(!std::isfinite(config.soft_penalty) || !(config.soft_penalty > scalar_type{0}))
        return ctrlpp::unexpected(moving_horizon_configuration_error::non_positive_soft_penalty);
    if(!std::isfinite(config.numerical_eps) || !(config.numerical_eps > scalar_type{0}))
        return ctrlpp::unexpected(moving_horizon_configuration_error::non_positive_numerical_eps);
    return {};
}

}

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N>
struct mhe_config
{
    using scalar_type = Scalar;

    static_assert(N > 0, "Window length N must be positive: it sizes the fixed estimation window arrays the estimator rotates and reads the trailing element of");

    Matrix<Scalar, NX, NX> Q{Matrix<Scalar, NX, NX>::Identity()};
    Matrix<Scalar, NY, NY> R{Matrix<Scalar, NY, NY>::Identity()};
    Vector<Scalar, NX> x0{Vector<Scalar, NX>::Zero()};
    Matrix<Scalar, NX, NX> P0{Matrix<Scalar, NX, NX>::Identity()};
    Scalar arrival_cost_weight{Scalar{1}};
    std::optional<Vector<Scalar, NX>> x_min{};
    std::optional<Vector<Scalar, NX>> x_max{};
    std::optional<Vector<Scalar, NY>> residual_bound{};
    bool soft_constraints{true};
    Scalar soft_penalty{Scalar{1e4}};
    Scalar numerical_eps{std::cbrt(std::numeric_limits<Scalar>::epsilon())};
};

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N, std::size_t NC = 0>
struct nmhe_config
{
    using scalar_type = Scalar;

    static_assert(N > 0, "Window length N must be positive: it sizes the fixed estimation window arrays the estimator rotates and reads the trailing element of");

    Matrix<Scalar, NX, NX> Q{Matrix<Scalar, NX, NX>::Identity()};
    Matrix<Scalar, NY, NY> R{Matrix<Scalar, NY, NY>::Identity()};
    Vector<Scalar, NX> x0{Vector<Scalar, NX>::Zero()};
    Matrix<Scalar, NX, NX> P0{Matrix<Scalar, NX, NX>::Identity()};
    Scalar arrival_cost_weight{Scalar{1}};
    std::optional<Vector<Scalar, NX>> x_min{};
    std::optional<Vector<Scalar, NX>> x_max{};
    std::optional<Vector<Scalar, NY>> residual_bound{};
    bool soft_constraints{true};
    Scalar soft_penalty{Scalar{1e4}};
    Scalar numerical_eps{std::cbrt(std::numeric_limits<Scalar>::epsilon())};
    std::optional<std::function<Vector<Scalar, NC>(const Vector<Scalar, NX>&)>> path_constraint_cb{};
    Vector<Scalar, NC> path_penalty{detail::default_penalty<Scalar, NC>()};
};

}

#endif
