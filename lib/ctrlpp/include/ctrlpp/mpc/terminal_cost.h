#ifndef HPP_GUARD_CTRLPP_MPC_TERMINAL_COST_H
#define HPP_GUARD_CTRLPP_MPC_TERMINAL_COST_H

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/control/dare.h"

#include "ctrlpp/mpc/qp_types.h"

#include "ctrlpp/model/state_space.h"

#include <Eigen/Dense>

#include <cstddef>
#include <utility>
#include <variant>
#include <optional>

namespace ctrlpp
{

/// @brief Construction error of the linear MPC factory: a configuration refusal, or, with Qf unset,
/// the Riccati solve's own refusal of the terminal cost it supplies.
using mpc_construction_error = std::variant<controller_construction_error, dare_error>;

}

namespace ctrlpp::detail
{

/// @brief Whether the Riccati solve would refuse R as singular.
///
/// This is the solve's own reciprocal-pivot test on its own operand, R divided by the weight scale it
/// equilibrates with. Non-finite weights are refused by the solve before that test, so they are no
/// verdict on R.
template <typename Scalar, std::size_t NX, std::size_t NU>
auto riccati_input_weight_singular(const Matrix<Scalar, NX, NU>& B, const Matrix<Scalar, NX, NX>& Q, const Matrix<Scalar, NU, NU>& R) -> bool
{
    if(!Q.allFinite() || !R.allFinite())
        return false;
    const Matrix<Scalar, NU, NU> R_scaled = R / dare_weight_scale<Scalar, NX, NU>(Q, R);
    auto const G = factor_dare_g<Scalar, NX, NU>(B, R_scaled);
    return !G.has_value() && G.error() == dare_error::singular_r;
}

template <typename Scalar, std::size_t NX, std::size_t NY>
struct mpc_terminal_weights
{
    Matrix<Scalar, NX, NX> state;
    Matrix<Scalar, NX, NY> linear;
};

/// @brief The linear MPC's terminal cost, in state space and as the output-space tracking weight.
///
/// With no terminal weight configured, the infinite-horizon Riccati solution is the terminal cost.
/// A singular R is refused first, because that solve inverts it, and a solve that declines for any
/// other cause is refused with its own enumerator: no other terminal cost is ever posed in its place.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
auto mpc_terminal_cost(const discrete_state_space<Scalar, NX, NU, NY>& system, const Matrix<Scalar, NX, NX>& Q_state, const Matrix<Scalar, NU, NU>& R,
                       const std::optional<Matrix<Scalar, NY, NY>>& Qf) -> expected<mpc_terminal_weights<Scalar, NX, NY>, mpc_construction_error>
{
    if(Qf.has_value())
        return mpc_terminal_weights<Scalar, NX, NY>{.state = system.C.transpose() * (*Qf) * system.C, .linear = system.C.transpose() * (*Qf)};
    if(riccati_input_weight_singular<Scalar, NX, NU>(system.B, Q_state, R))
        return ctrlpp::unexpected(mpc_construction_error{controller_construction_error::singular_input_weight});

    auto riccati = dare<Scalar, NX, NU>(system.A, system.B, Q_state, R);
    if(!riccati)
        return ctrlpp::unexpected(mpc_construction_error{riccati.error()});
    // The tracking term maps the state weight to output space through the pseudo-inverse C^T (C C^T)^{-1}.
    Matrix<Scalar, NY, NY> CCt = system.C * system.C.transpose();
    Matrix<Scalar, NX, NY> linear = riccati->P * system.C.transpose() * CCt.colPivHouseholderQr().solve(Matrix<Scalar, NY, NY>::Identity());
    return mpc_terminal_weights<Scalar, NX, NY>{.state = std::move(riccati->P), .linear = std::move(linear)};
}

}

#endif
