#ifndef HPP_GUARD_CTRLPP_MPC_H
#define HPP_GUARD_CTRLPP_MPC_H

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/qp_solver.h"
#include "ctrlpp/mpc/mpc_config.h"
#include "ctrlpp/mpc/diagnostics.h"
#include "ctrlpp/mpc/terminal_set.h"
#include "ctrlpp/mpc/terminal_cost.h"
#include "ctrlpp/mpc/qp_formulation.h"
#include "ctrlpp/mpc/config_validation.h"

#include "ctrlpp/mpc/detail/mpc_solve.h"

#include "ctrlpp/model/state_space.h"

#include <Eigen/Dense>

#include <span>
#include <limits>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <variant>
#include <type_traits>

namespace ctrlpp
{

template <typename Scalar, std::size_t NX, std::size_t NU, qp_solver Solver, std::size_t NY = NX>
class mpc : private detail::mpc_solve<Scalar, NX, NU, Solver, NY, mpc_config<Scalar, NX, NU, NY>>
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);

    using setup_error = detail::qp_setup_error_t<Solver, qp_problem<Scalar>>;
    using base = detail::mpc_solve<Scalar, NX, NU, Solver, NY, mpc_config<Scalar, NX, NU, NY>>;

    // Passkey: only the controller can make one, so create stays the only way to originate a controller and its
    // validation and solver setup cannot be bypassed.
    class validated_tag
    {
        friend mpc;
        validated_tag() = default;
    };

public:
    using construction_error = mpc_construction_error<setup_error>;

    // Equal alternative types would make the variant's construction by type ambiguous.
    static_assert(!std::is_same_v<setup_error, controller_construction_error> && !std::is_same_v<setup_error, dare_error>,
            "the solver's setup-error type must differ from controller_construction_error and dare_error");

    // Public so ctrlpp::expected can construct the controller in the caller's
    // storage; unreachable without the validated_tag only the controller can make.
    mpc(validated_tag, const discrete_state_space<Scalar, NX, NU, NY>& system, const mpc_config<Scalar, NX, NU, NY>& config, Matrix<Scalar, NX, NX> Q_state,
            detail::mpc_terminal_weights<Scalar, NX, NY> terminal, Solver&& solver)
        : base{system, config, std::move(Q_state), std::move(terminal), std::move(solver)}
    {
    }

    static expected<mpc, construction_error> create(const discrete_state_space<Scalar, NX, NU, NY>& system, const mpc_config<Scalar, NX, NU, NY>& config)
    {
        return create(system, config, Solver{});
    }

    // The solver is moved in before the initial QP is posed, so its settings govern setup.
    static expected<mpc, construction_error> create(const discrete_state_space<Scalar, NX, NU, NY>& system, const mpc_config<Scalar, NX, NU, NY>& config, Solver solver)
    {
        if(config.horizon <= 0)
            return unexpected(construction_error{controller_construction_error::non_positive_horizon});
        if(config.horizon > horizon_bound(config))
            return unexpected(construction_error{controller_construction_error::horizon_overflow});
        auto admitted = detail::admit_mpc_config(config);
        if(!admitted)
            return unexpected(construction_error{admitted.error()});
        if(!system.A.allFinite() || !system.B.allFinite() || !system.C.allFinite())
            return unexpected(construction_error{controller_construction_error::non_finite_system});
        Matrix<Scalar, NX, NX> Q_state = system.C.transpose() * admitted->Q * system.C;
        if(!Q_state.allFinite())
            return unexpected(construction_error{controller_construction_error::non_finite_state_space_weight});
        auto terminal = detail::mpc_terminal_cost<setup_error, Scalar, NX, NU, NY>(system, Q_state, admitted->R, admitted->Qf);
        if(!terminal)
            return unexpected(terminal.error());
        return construct_and_set_up(system, *admitted, std::move(Q_state), std::move(*terminal), std::move(solver));
    }

    expected<solve_output<Scalar, NU>, solver_error> solve(const Vector<Scalar, NX>& x0)
    {
        this->m_update.q.setZero();
        return this->solve_impl(x0);
    }

    expected<solve_output<Scalar, NU>, solver_error> solve(const Vector<Scalar, NX>& x0, const Vector<Scalar, NY>& y_ref)
    {
        std::int32_t N = this->m_config.horizon;
        Vector<Scalar, NX> CtQ_yref = this->m_CtQ * y_ref;
        for(std::int32_t k = 0; k < N; ++k)
            this->m_update.q.segment(k * nx, nx) = -CtQ_yref;
        this->m_update.q.segment(N * nx, nx) = -(this->m_Qf_linear * y_ref);
        this->m_update.q.segment(this->m_n_x_total, this->m_n_dec - this->m_n_x_total).setZero();
        return this->solve_impl(x0);
    }

    expected<solve_output<Scalar, NU>, solver_error> solve(const Vector<Scalar, NX>& x0, std::span<const Vector<Scalar, NY>> y_ref)
    {
        std::int32_t N = this->m_config.horizon;
        if(y_ref.size() < static_cast<std::size_t>(N + 1))
            return unexpected(solver_error::invalid_problem);
        for(std::int32_t k = 0; k < N; ++k)
            this->m_update.q.segment(k * nx, nx) = -(this->m_CtQ * y_ref[static_cast<std::size_t>(k)]);
        this->m_update.q.segment(N * nx, nx) = -(this->m_Qf_linear * y_ref[static_cast<std::size_t>(N)]);
        this->m_update.q.segment(this->m_n_x_total, this->m_n_dec - this->m_n_x_total).setZero();
        return this->solve_impl(x0);
    }

    void set_applied_input(const Vector<Scalar, NU>& u) { this->m_u_prev = u; }

    expected<std::pair<std::vector<Vector<Scalar, NX>>, std::vector<Vector<Scalar, NU>>>, solver_error> trajectory() const
    {
        if(!this->m_has_solution)
            return unexpected(solver_error::setup_incomplete);

        // No length check: a solution is stored only after extract_solution has checked that the primal covers the
        // decision dimension these offsets come from.
        std::int32_t N = this->m_config.horizon;
        std::vector<Vector<Scalar, NX>> states;
        std::vector<Vector<Scalar, NU>> inputs;
        states.reserve(static_cast<std::size_t>(N + 1));
        inputs.reserve(static_cast<std::size_t>(N));

        for(std::int32_t k = 0; k <= N; ++k)
            states.push_back(this->m_last_primal.segment(k * nx, nx));
        for(std::int32_t k = 0; k < N; ++k)
            inputs.push_back(this->m_last_primal.segment(this->m_n_x_total + k * nu, nu));

        return std::pair{std::move(states), std::move(inputs)};
    }

    mpc_diagnostics<Scalar> diagnostics() const { return this->m_last_diagnostics; }

private:
    // The setup failure is assigned into the one named result, and that result is
    // this function's only return: gcc elides a named result only when every
    // return names it, so the configuration refusals live in create instead.
    static expected<mpc, construction_error> construct_and_set_up(const discrete_state_space<Scalar, NX, NU, NY>& system, const mpc_config<Scalar, NX, NU, NY>& config,
            Matrix<Scalar, NX, NX>&& Q_state, detail::mpc_terminal_weights<Scalar, NX, NY>&& terminal, Solver&& solver)
    {
        expected<mpc, construction_error> result{std::in_place, validated_tag{}, system, config, std::move(Q_state), std::move(terminal), std::move(solver)};
        auto setup = detail::setup_qp_solver(result->m_solver, result->initial_problem());
        if(!setup.has_value())
            result = unexpected(construction_error{setup.error()});
        return result;
    }

    // Both derived dimensions, N (2 nx + nu) + nx decision variables and N (2 nx + 2 nu) + nx + n_terminal constraint
    // rows at the worst-case configuration, stay representable in the horizon's own type up to this horizon.
    static std::int32_t horizon_bound(const mpc_config<Scalar, NX, NU, NY>& config)
    {
        constexpr std::int32_t per_step = 2 * nx + 2 * nu;
        const std::int32_t constant_dimensions = nx + static_cast<std::int32_t>(detail::terminal_constraint_rows<Scalar, NX>(config.terminal_constraint_set));
        return (std::numeric_limits<std::int32_t>::max() - constant_dimensions) / per_step;
    }
};

}

#endif
