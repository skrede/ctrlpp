#ifndef HPP_GUARD_CTRLPP_NMPC_H
#define HPP_GUARD_CTRLPP_NMPC_H

#include "ctrlpp/expected.h"

#include "ctrlpp/model/dynamics_model.h"

#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/nmpc_config.h"
#include "ctrlpp/mpc/nlp_formulation.h"
#include "ctrlpp/mpc/config_validation.h"

#include "ctrlpp/mpc/detail/nmpc_static_solve.h"
#include "ctrlpp/mpc/detail/nmpc_dynamic_solve.h"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <type_traits>

namespace ctrlpp
{

// The horizon is known only at construction, so the solve binds a runtime-erased nlp_problem and the solver's
// allocating path. The allocation-free default is nmpc_static.
template <typename Scalar, std::size_t NX, std::size_t NU, nlp_solver Solver, dynamics_model<Scalar, NX, NU> Dynamics, std::size_t NC = 0, std::size_t NTC = 0>
class nmpc_dynamic : private detail::nmpc_dynamic_solve<Scalar, NX, NU, Solver, Dynamics, NC, NTC>
{
    using base = detail::nmpc_dynamic_solve<Scalar, NX, NU, Solver, Dynamics, NC, NTC>;
    using setup_error = detail::nlp_setup_error_t<Solver, nlp_problem<Scalar>>;

    // Only the controller can make the passkey, so create stays the only way to originate a controller and the
    // validation and the solver setup cannot be bypassed.
    class validated_tag
    {
        friend nmpc_dynamic;
        validated_tag() = default;
    };

public:
    using construction_error = nmpc_construction_error<setup_error>;

    // Equal alternative types would make the variant's construction by type ambiguous.
    static_assert(!std::is_same_v<setup_error, controller_construction_error> && !std::is_same_v<setup_error, nlp_formulation_error>,
            "the solver's setup-error type must differ from controller_construction_error and nlp_formulation_error");

    // Public so ctrlpp::expected can construct the controller in the caller's storage; unreachable without the
    // passkey.
    nmpc_dynamic(validated_tag, Dynamics dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config, Solver&& solver)
        : base{std::move(dynamics), config, std::move(solver)}
    {
    }

    static expected<nmpc_dynamic, construction_error> create(Dynamics dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config)
    {
        return create(std::move(dynamics), config, Solver{});
    }

    static expected<nmpc_dynamic, construction_error> create(Dynamics dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config, Solver solver)
    {
        if(config.horizon <= 0)
            return unexpected(construction_error{controller_construction_error::non_positive_horizon});
        if(config.horizon > base::horizon_bound())
            return unexpected(construction_error{controller_construction_error::horizon_overflow});
        auto admitted = detail::admit_nmpc_config(config);
        if(!admitted)
            return unexpected(construction_error{admitted.error()});
        return construct_and_set_up(std::move(dynamics), *admitted, std::move(solver));
    }

    using base::solve;
    using base::set_applied_input;
    using base::trajectory;
    using base::diagnostics;
    using base::problem;
    using base::last_solution;

private:
    // The setup failure is assigned into the one named result, its only return: gcc elides a named result only
    // when every return names it, so the refusals live in create instead.
    static expected<nmpc_dynamic, construction_error> construct_and_set_up(Dynamics&& dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config, Solver&& solver)
    {
        expected<nmpc_dynamic, construction_error> result{std::in_place, validated_tag{}, std::move(dynamics), config, std::move(solver)};
        auto setup = detail::setup_nlp_solver(result->m_solver, *result->m_problem);
        if(!setup.has_value())
            result = unexpected(construction_error{setup.error()});
        return result;
    }
};

// The horizon NH is a template parameter, so the decision dimension NV = (NH + 1) NX + NH NU is a compile-time
// constant and the steady-state solve on argmin's fixed-N floor allocates nothing. Only the slack-free formulation is
// posed: create refuses a horizon other than NH and a configuration that would add slack.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NH, typename Solver, dynamics_model<Scalar, NX, NU> Dynamics, std::size_t NC = 0, std::size_t NTC = 0>
class nmpc_static : private detail::nmpc_static_solve<Scalar, NX, NU, NH, Solver, Dynamics, NC, NTC>
{
    static_assert(NH > 0,
            "Horizon NH must be positive: a zero horizon leaves no input to apply and makes the horizon - 1 warm-start shift and the input offset into the primal read outside the "
            "decision vector");

    using base = detail::nmpc_static_solve<Scalar, NX, NU, NH, Solver, Dynamics, NC, NTC>;
    using config_type = nmpc_config<Scalar, NX, NU, NC, NTC>;

public:
    using base::horizon;
    using base::problem_dimension;
    using typename base::problem_type;

private:
    using setup_error = detail::nlp_setup_error_t<Solver, problem_type>;

    // Only the controller can make the passkey, so create stays the only way to originate a controller and the
    // validation, the formulation refusals and the solver setup cannot be bypassed.
    class validated_tag
    {
        friend nmpc_static;
        validated_tag() = default;
    };

public:
    using construction_error = nmpc_construction_error<setup_error>;

    // Equal alternative types would make the variant's construction by type ambiguous.
    static_assert(!std::is_same_v<setup_error, controller_construction_error> && !std::is_same_v<setup_error, nlp_formulation_error>,
            "the solver's setup-error type must differ from controller_construction_error and nlp_formulation_error");

    // Public so ctrlpp::expected can construct the controller in the caller's storage; unreachable without the
    // passkey. The solver is built in place: a solver that carries its workspace inline (argmin's fixed-N policies,
    // tens of kilobytes) would otherwise exist twice.
    nmpc_static(validated_tag, Dynamics dynamics, const config_type& config)
        : base{std::move(dynamics), config}
    {
    }

    // Taken by rvalue reference: a by-value parameter would hold a second solver in the in-place construction's own
    // frame at -O0.
    nmpc_static(validated_tag, Dynamics dynamics, const config_type& config, Solver&& solver)
        : base{std::move(dynamics), config, std::move(solver)}
    {
    }

    static expected<nmpc_static, construction_error> create(Dynamics dynamics, const config_type& config) { return admit_and_construct(std::move(dynamics), config); }

    static expected<nmpc_static, construction_error> create(Dynamics dynamics, const config_type& config, Solver solver)
    {
        return admit_and_construct(std::move(dynamics), config, std::move(solver));
    }

    using base::solve;
    using base::set_applied_input;
    using base::diagnostics;
    using base::problem;
    using base::last_solution;

private:
    // The pack is empty for the default solver and one solver otherwise, so both overloads share one admission
    // sequence without a default-built solver.
    template <typename... SolverArgument>
    static expected<nmpc_static, construction_error> admit_and_construct(Dynamics&& dynamics, const config_type& config, SolverArgument&&... solver)
    {
        auto admitted = detail::admit_nmpc_config(config);
        if(!admitted.has_value())
            return unexpected(construction_error{admitted.error()});
        if(auto refused = detail::static_formulation_refusal<NH>(*admitted); !refused.has_value())
            return unexpected(construction_error{refused.error()});
        return construct_and_set_up(std::move(dynamics), *admitted, std::forward<SolverArgument>(solver)...);
    }

    // The posing and setup refusals are assigned into the one named result, its only return: gcc elides a named
    // result only when every return names it. Posing repeats the formulation refusals the caller already ran, so a
    // divergence between the two is refused here, not set up against an empty problem.
    template <typename... SolverArgument>
    static expected<nmpc_static, construction_error> construct_and_set_up(Dynamics&& dynamics, const config_type& config, SolverArgument&&... solver)
    {
        expected<nmpc_static, construction_error> result{std::in_place, validated_tag{}, std::move(dynamics), config, std::forward<SolverArgument>(solver)...};
        if(auto posed = result->pose_problem(); !posed.has_value())
            result = unexpected(construction_error{posed.error()});
        else if(auto setup = detail::setup_nlp_solver(result->m_solver, *result->m_problem); !setup.has_value())
            result = unexpected(construction_error{setup.error()});
        return result;
    }
};

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NH, typename Solver, dynamics_model<Scalar, NX, NU> Dynamics, std::size_t NC = 0, std::size_t NTC = 0>
using nmpc = nmpc_static<Scalar, NX, NU, NH, Solver, Dynamics, NC, NTC>;

}

#endif
