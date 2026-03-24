#ifndef HPP_GUARD_CTRLPP_MHE_MHE_CONFIG_H
#define HPP_GUARD_CTRLPP_MHE_MHE_CONFIG_H

#include "ctrlpp/mpc/nmpc_config.h"
#include "ctrlpp/types.h"

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <optional>

namespace ctrlpp
{

/// Configuration for linear MHE (QP-based).
///
/// Template parameters:
///   Scalar - floating-point type
///   NX     - state dimension
///   NU     - input dimension
///   NY     - measurement dimension
///   N      - horizon length (number of measurement steps in window)
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N>
struct mhe_config
{
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
    Scalar numerical_eps{std::sqrt(std::numeric_limits<Scalar>::epsilon())};
};

/// Configuration for nonlinear MHE (NLP-based).
///
/// Extends mhe_config with general nonlinear path constraints g(x) <= 0.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N, std::size_t NC = 0>
struct nmhe_config
{
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
    Scalar numerical_eps{std::sqrt(std::numeric_limits<Scalar>::epsilon())};
    std::optional<std::function<Vector<Scalar, NC>(const Vector<Scalar, NX>&)>> path_constraint{};
    Vector<Scalar, NC> path_penalty{detail::default_penalty<Scalar, NC>()};
};

} // namespace ctrlpp

#endif
