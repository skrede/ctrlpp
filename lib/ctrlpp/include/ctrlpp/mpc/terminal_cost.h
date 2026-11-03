#ifndef HPP_GUARD_CTRLPP_MPC_TERMINAL_COST_H
#define HPP_GUARD_CTRLPP_MPC_TERMINAL_COST_H

#include "ctrlpp/mpc/qp_types.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/control/dare.h"

#include "ctrlpp/model/state_space.h"

#include <Eigen/Dense>

#include <cstddef>
#include <utility>
#include <variant>
#include <optional>

namespace ctrlpp
{

// The dare_error alternative arises only with Qf unset; SetupError is the bound solver's own setup error.
template <typename SetupError>
using mpc_construction_error = std::variant<controller_construction_error, dare_error, SetupError>;

}

namespace ctrlpp::detail
{

// The Riccati solve's own reciprocal-pivot test on R divided by its equilibration scale. Non-finite
// weights are no verdict on R: the solve refuses them before that test.
template <typename Scalar, std::size_t NX, std::size_t NU>
bool riccati_input_weight_singular(const Matrix<Scalar, NX, NU>& B, const Matrix<Scalar, NX, NX>& Q, const Matrix<Scalar, NU, NU>& R)
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

// With Qf unset the infinite-horizon Riccati solution is the terminal cost, and a singular R is refused first
// because that solve inverts it. A declined solve is refused with its own enumerator, never replaced.
template <typename SetupError, typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
expected<mpc_terminal_weights<Scalar, NX, NY>, mpc_construction_error<SetupError>> mpc_terminal_cost(
        const discrete_state_space<Scalar, NX, NU, NY>& system, const Matrix<Scalar, NX, NX>& Q_state, const Matrix<Scalar, NU, NU>& R, const std::optional<Matrix<Scalar, NY, NY>>& Qf)
{
    using construction_error = mpc_construction_error<SetupError>;
    if(Qf.has_value())
    {
        Matrix<Scalar, NX, NX> state = system.C.transpose() * (*Qf) * system.C;
        if(!state.allFinite())
            return ctrlpp::unexpected(construction_error{controller_construction_error::non_finite_state_space_weight});
        return mpc_terminal_weights<Scalar, NX, NY>{.state = std::move(state), .linear = system.C.transpose() * (*Qf)};
    }
    if(riccati_input_weight_singular<Scalar, NX, NU>(system.B, Q_state, R))
        return ctrlpp::unexpected(construction_error{controller_construction_error::singular_input_weight});

    auto riccati = dare<Scalar, NX, NU>(system.A, system.B, Q_state, R);
    if(!riccati)
        return ctrlpp::unexpected(construction_error{riccati.error()});
    // The tracking term maps the state weight to output space through the pseudo-inverse C^T (C C^T)^{-1}.
    Matrix<Scalar, NY, NY> CCt = system.C * system.C.transpose();
    Matrix<Scalar, NX, NY> linear = riccati->P * system.C.transpose() * CCt.colPivHouseholderQr().solve(Matrix<Scalar, NY, NY>::Identity());
    return mpc_terminal_weights<Scalar, NX, NY>{.state = std::move(riccati->P), .linear = std::move(linear)};
}

}

#endif
