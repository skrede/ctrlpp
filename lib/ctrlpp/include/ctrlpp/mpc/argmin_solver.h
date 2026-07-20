#ifndef HPP_GUARD_CTRLPP_MPC_ARGMIN_SOLVER_H
#define HPP_GUARD_CTRLPP_MPC_ARGMIN_SOLVER_H

#include "ctrlpp/config.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/mpc/nlp_types.h"
#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/argmin_problem.h"
#include "ctrlpp/mpc/argmin_policies.h"

#include <argmin/result/status.h>

#include <argmin/solver/options.h>
#include <argmin/solver/step_budget_solver.h>
#include <argmin/solver/time_budget_options.h>
#include <argmin/solver/step_and_time_budget_solver.h>

#include <chrono>
#include <memory>
#include <cstdint>
#include <variant>
#include <stdexcept>
#include <type_traits>

namespace ctrlpp
{

// argmin_solver<Scalar, Policy, Constrained, NV>: the trailing NV selects the
// solve path. NV == Eigen::Dynamic (the DEFAULT) is the original runtime-erased
// solver, byte-for-byte: it binds an nlp_problem<Scalar>, instantiates argmin's
// dynamic step_budget_solver, and every existing caller compiles unchanged.
// NV != Eigen::Dynamic is the allocation-free static path (Route A / SEED-002):
// the policy algorithm is rebound to the compile-time dimension NV, argmin's
// compile-time-N step_budget_solver is instantiated, the bridge uses fixed-size
// decision-vector storage, and setup binds an nlp_problem_static<Scalar, NV>.
template <typename Scalar, typename Policy, bool Constrained = true, int NV = Eigen::Dynamic>
class argmin_solver
{
    // Compile-time guard for the structurally-always-invalid configuration: a
    // raw MMA-family policy (argmin_mma / argmin_gcmma) driving the constrained
    // bridge. NMPC always emits equality constraints (initial-state pin +
    // continuity), so a constrained raw-MMA solver can never solve the problem
    // it is handed. This is the compile-time complement to try_setup's runtime
    // reject: the static_assert bars the always-wrong constrained NMPC
    // instantiation, while try_setup rejects a general nlp_problem that happens
    // to carry equalities. Wrap the policy as argmin_auglag<argmin_mma> (the
    // outer augmented-Lagrangian loop absorbs the equalities) to fix it.
    static_assert(!(is_raw_mma_family_v<Policy> && Constrained),
        "Raw MMA-family policies (argmin_mma / argmin_gcmma) cannot handle "
        "equality constraints and NMPC always emits them; use "
        "argmin_auglag<argmin_mma> for constrained problems, or set "
        "Constrained=false for unconstrained (bound-only) use.");

public:
    using scalar_type = Scalar;
    // On the static path (NV != Eigen::Dynamic) the policy algorithm is rebound
    // to the compile-time dimension so argmin's fixed-N NW-SQP state sizes its
    // decision-vector buffers with fixed-size Eigen types. On the dynamic default
    // the rebind resolves to the same type as Policy::algorithm (byte-identical).
    using argmin_policy = rebind_argmin_algorithm_t<typename Policy::algorithm, NV>;
    using bridge_type = std::conditional_t<Constrained,
        argmin_constrained_problem<Scalar, NV>,
        argmin_problem<Scalar, NV>>;
    // The bound problem type follows NV: nlp_problem_static on the static path,
    // the runtime-erased nlp_problem on the dynamic default.
    using problem_type = typename bridge_type::problem_type;
    using step_solver_type =
        argmin::step_budget_solver<argmin_policy, NV, bridge_type,
            argmin_ctrlpp_convergence>;
    using timed_solver_type =
        argmin::step_and_time_budget_solver<argmin_policy, NV, bridge_type,
            argmin_ctrlpp_convergence>;
    using solver_storage_type =
        std::variant<std::monostate, step_solver_type, timed_solver_type>;
    using settings_type = std::conditional_t<
        is_mma_family_v<Policy>,
        argmin_mma_settings<Scalar>,
        argmin_settings<Scalar>>;

    explicit argmin_solver(settings_type settings = {})
        : settings_{settings}
        , bridge_{std::make_unique<bridge_type>()}
    {}

    // Move is correct-by-default: `bridge_` lives behind a `unique_ptr` (a stable
    // heap address), so moving relocates only the owning pointer while the
    // pointed-to bridge stays put. argmin's `solver_core` caches the problem BY
    // REFERENCE (`const Problem* problem_ptr_`, solver_core.h:657) and its policy
    // state likewise caches `&problem`; both back-pointers therefore remain valid
    // across a defaulted move. A hand-written move ctor is deliberately avoided.
    argmin_solver(argmin_solver&&) = default;
    argmin_solver& operator=(argmin_solver&&) = default;

    // Fork-copy: the ONE deliberately hand-written special member. argmin's
    // `solver_core` is copy-deleted (solver_core.h:306-307), so a defaulted copy
    // is ill-formed. This snapshots the settings and the external problem pointer,
    // allocates a FRESH bridge rebound to the same problem, and leaves `solver_`
    // as `std::monostate` so it is re-emplaced lazily on the next solve()/step()
    // via prepare_solver's monostate branch. The result is a clean independent
    // fork, matching the intended clean-snapshot copy semantics, with no argmin copy support.
    argmin_solver(const argmin_solver& other)
        : settings_{other.settings_}
        , problem_{other.problem_}
        , bridge_{std::make_unique<bridge_type>()}
        , solver_{}
    {
        rebind_bridge();
    }

    argmin_solver& operator=(const argmin_solver& other)
    {
        if(this != &other)
        {
            auto fresh = std::make_unique<bridge_type>();
            settings_ = other.settings_;
            problem_ = other.problem_;
            // Drop the stale solver (which caches the old bridge address) BEFORE
            // swapping in the fresh bridge, so no argmin back-pointer ever
            // observes a freed object.
            solver_.template emplace<std::monostate>();
            bridge_ = std::move(fresh);
            rebind_bridge();
        }
        return *this;
    }

    /// @brief Fallible setup: binds the problem into the bridge and rejects a
    /// raw MMA-family policy handed equality constraints. Returns an empty
    /// expected on success and `argmin_setup_error::incompatible_equality_constraints`
    /// when a raw MMA/GCMMA policy is given a problem that carries equalities.
    /// Reported through the `ctrlpp::expected` channel so the reject works in
    /// all build modes, including `-fno-exceptions`.
    [[nodiscard]] auto try_setup(const problem_type& problem)
        -> ctrlpp::expected<void, argmin_setup_error>
    {
        problem_ = &problem;
        solver_.template emplace<std::monostate>();

        if constexpr(Constrained)
            bridge_->partition(problem);
        else
            bridge_->bind(problem);

        if constexpr(is_raw_mma_family_v<Policy>)
        {
            // Raw MMA/GCMMA cannot represent equality constraints (the
            // auglag-wrapped variants absorb them and are exempt via
            // is_raw_mma_family_v). Scan the raw problem for equalities the
            // same way nlopt_solver does, rather than the constrained bridge's
            // num_equality(): the class-body static_assert bars the constrained
            // raw-MMA bridge, so the only compilable raw-MMA configuration is
            // Constrained=false, whose bridge does not partition constraints.
            // Rejecting an equality-carrying problem here prevents a
            // silently-wrong solve that would ignore those equalities.
            if(problem_has_equality(problem))
                return ctrlpp::unexpected(argmin_setup_error::incompatible_equality_constraints);
        }

        return {};
    }

#if CTRLPP_HAS_EXCEPTIONS
    /// @brief Throwing convenience wrapper over `try_setup`. Throws
    /// `std::invalid_argument` when a raw MMA-family policy is handed equality
    /// constraints.
    void setup(const problem_type& problem)
    {
        if(try_setup(problem).has_value())
            return;

        throw std::invalid_argument(
            "raw MMA-family argmin policy cannot handle equality constraints; "
            "wrap it as argmin_auglag<argmin_mma>");
    }
#endif

    auto solve(const nlp_update<Scalar>& update) -> nlp_result<Scalar>
    {
        prepare_solver(update.x0);
        return with_solver([&](auto& solver)
        {
            auto result = solver.solve();
            return translate_result(result);
        });
    }

    auto step(const nlp_update<Scalar>& update, int max_steps) -> nlp_result<Scalar>
    {
        prepare_solver(update.x0);
        return with_solver([&](auto& solver)
        {
            auto result = solver.step_n(static_cast<std::uint32_t>(max_steps));
            return translate_result(result);
        });
    }

private:
    // Rebind the (fresh) bridge to the retained external problem. Shared by the
    // fork-copy ctor and copy-assignment. A null problem_ (never set up) leaves
    // the bridge unbound, exactly as a default-constructed solver would be.
    void rebind_bridge()
    {
        if(problem_ == nullptr)
            return;

        if constexpr(Constrained)
            bridge_->partition(*problem_);
        else
            bridge_->bind(*problem_);
    }

    static auto problem_has_equality(const problem_type& problem) -> bool
    {
        for(int i = 0; i < problem.n_constraints; ++i)
        {
            if(problem.c_lower[i] == problem.c_upper[i])
                return true;
        }
        return false;
    }

    void prepare_solver(const Eigen::VectorX<Scalar>& x0)
    {
        if(std::holds_alternative<std::monostate>(solver_))
        {
            if(has_time_budget())
                emplace_timed_solver(x0);
            else
                emplace_step_solver(x0);

            return;
        }

        const auto ws = [&]() -> warm_start_mode
        {
            if constexpr(is_mma_family_v<Policy>)
                return settings_.base.warm_start;
            else
                return settings_.warm_start;
        }();

        std::visit([&](auto& solver)
        {
            using solver_t = std::decay_t<decltype(solver)>;
            if constexpr(!std::is_same_v<solver_t, std::monostate>)
            {
                if(ws == warm_start_mode::curvature)
                    solver.reset(x0);
                else
                    solver.reset_clear(x0);
            }
        }, solver_);
    }

    auto make_solver_options() const -> argmin::solver_options<argmin_ctrlpp_convergence>
    {
        argmin::solver_options<argmin_ctrlpp_convergence> opts;

        auto const& s = [&]() -> auto const&
        {
            if constexpr(is_mma_family_v<Policy>)
                return settings_.base;
            else
                return settings_;
        }();

        opts.max_iterations = static_cast<std::uint32_t>(s.max_eval);

        if constexpr(Constrained)
        {
            if(s.constraint_tol > Scalar{0})
            {
                opts.constraint_tolerance = static_cast<double>(s.constraint_tol);
                opts.feasibility_tolerance = static_cast<double>(s.constraint_tol);
            }
        }

        // ftol_rel / xtol_rel are relative tolerances: wire them to argmin's
        // RELATIVE criteria via the _rel setters. argmin_ctrlpp_convergence
        // carries objective_tolerance_rel_criterion / step_tolerance_rel_criterion
        // so these requires-guarded setters are well-formed. This is a
        // deliberate convergence-behavior change from argmin's default absolute
        // criteria (see argmin_ctrlpp_convergence in argmin_policies.h).
        opts.set_objective_threshold_rel(static_cast<double>(s.ftol_rel));
        opts.set_step_threshold_rel(static_cast<double>(s.xtol_rel));

        return opts;
    }

    auto make_time_budget_options() const -> argmin::time_budget_options<argmin_ctrlpp_convergence>
    {
        argmin::time_budget_options<argmin_ctrlpp_convergence> opts;
        opts.core = make_solver_options();

        auto const& s = [&]() -> auto const&
        {
            if constexpr(is_mma_family_v<Policy>)
                return settings_.base;
            else
                return settings_;
        }();

        opts.max_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::duration<double>(static_cast<double>(s.max_time)));

        return opts;
    }

    auto has_time_budget() const -> bool
    {
        auto const& s = [&]() -> auto const&
        {
            if constexpr(is_mma_family_v<Policy>)
                return settings_.base;
            else
                return settings_;
        }();

        return s.max_time > Scalar{0};
    }

    void emplace_step_solver(const Eigen::VectorX<Scalar>& x0)
    {
        if constexpr(is_mma_family_v<Policy>)
        {
            solver_.template emplace<step_solver_type>(
                argmin_policy{}, *bridge_, x0, make_solver_options(), make_mma_policy_opts());
        }
        else
        {
            solver_.template emplace<step_solver_type>(
                argmin_policy{}, *bridge_, x0, make_solver_options());
        }
    }

    void emplace_timed_solver(const Eigen::VectorX<Scalar>& x0)
    {
        if constexpr(is_mma_family_v<Policy>)
        {
            solver_.template emplace<timed_solver_type>(
                argmin_policy{}, *bridge_, x0, make_time_budget_options(), make_mma_policy_opts());
        }
        else
        {
            solver_.template emplace<timed_solver_type>(
                argmin_policy{}, *bridge_, x0, make_time_budget_options());
        }
    }

    auto make_mma_policy_opts() const -> typename argmin_policy::options_type
    {
        typename argmin_policy::options_type po{};
        if constexpr(std::is_same_v<Policy, argmin_mma>)
        {
            po.asymptote_init = static_cast<double>(settings_.asymptote_init);
            po.asymptote_expand = static_cast<double>(settings_.asymptote_incr);
            po.asymptote_contract = static_cast<double>(settings_.asymptote_decr);
        }
        else if constexpr(std::is_same_v<Policy, argmin_gcmma>)
        {
            po.asymptote_init = static_cast<double>(settings_.asymptote_init);
            po.asymptote_expand = static_cast<double>(settings_.asymptote_incr);
            po.asymptote_contract = static_cast<double>(settings_.asymptote_decr);
            po.max_inner_iterations =
                static_cast<std::uint16_t>(settings_.gcmma_inner_max);
        }
        else
        {
            // argmin_auglag<Inner> with is_mma_family_v<Inner>.
            if constexpr(std::is_same_v<Policy, argmin_auglag<argmin_mma>>)
            {
                po.inner_opts.asymptote_init =
                    static_cast<double>(settings_.asymptote_init);
                po.inner_opts.asymptote_expand =
                    static_cast<double>(settings_.asymptote_incr);
                po.inner_opts.asymptote_contract =
                    static_cast<double>(settings_.asymptote_decr);
            }
            else if constexpr(std::is_same_v<Policy, argmin_auglag<argmin_gcmma>>)
            {
                po.inner_opts.asymptote_init =
                    static_cast<double>(settings_.asymptote_init);
                po.inner_opts.asymptote_expand =
                    static_cast<double>(settings_.asymptote_incr);
                po.inner_opts.asymptote_contract =
                    static_cast<double>(settings_.asymptote_decr);
                po.inner_opts.max_inner_iterations =
                    static_cast<std::uint16_t>(settings_.gcmma_inner_max);
            }
        }
        return po;
    }

    static constexpr auto map_status(argmin::solver_status s) -> solve_status
    {
        switch(s)
        {
        case argmin::solver_status::converged:
        case argmin::solver_status::ftol_reached:
        case argmin::solver_status::xtol_reached:
            return solve_status::optimal;
        case argmin::solver_status::max_iterations:
        case argmin::solver_status::budget_exhausted:
        case argmin::solver_status::maxeval_reached:
            return solve_status::max_iterations;
        case argmin::solver_status::time_limit_reached:
            return solve_status::time_limit;
        case argmin::solver_status::stalled:
        case argmin::solver_status::roundoff_limited:
        case argmin::solver_status::objective_stalled:
        case argmin::solver_status::trust_region_step_rejected:
            return solve_status::solved_inaccurate;
        case argmin::solver_status::invalid_problem:
        case argmin::solver_status::diverged:
        case argmin::solver_status::aborted:
        case argmin::solver_status::running:
            return solve_status::error;
        }
        return solve_status::error;
    }

    template <typename Result>
    auto translate_result(const Result& r) const -> nlp_result<Scalar>
    {
        return nlp_result<Scalar>{
            .status = map_status(r.status),
            .x = r.x,
            .objective = r.objective_value,
            .solve_time = result_wall_time(r),
            .iterations = static_cast<int>(r.iterations),
            .primal_residual = r.constraint_violation,
        };
    }

    template <typename Result>
    static auto result_wall_time(const Result& r) -> Scalar
    {
        if constexpr(requires { r.wall_time; })
        {
            return static_cast<Scalar>(std::chrono::duration<double>(r.wall_time).count());
        }
        else
        {
            return Scalar{0};
        }
    }

    template <typename SolverFn>
    auto with_solver(SolverFn&& fn) -> nlp_result<Scalar>
    {
        return std::visit([&](auto& solver) -> nlp_result<Scalar>
        {
            using solver_t = std::decay_t<decltype(solver)>;
            if constexpr(std::is_same_v<solver_t, std::monostate>)
            {
                return nlp_result<Scalar>{.status = solve_status::error};
            }
            else
            {
                return fn(solver);
            }
        }, solver_);
    }

    settings_type settings_;
    const problem_type* problem_{nullptr};
    // Held behind a stable heap address so a defaulted move keeps argmin's
    // by-reference problem back-pointer (solver_core.h:657) valid. Declared
    // before `solver_` so that, at destruction, `solver_` (which caches this
    // bridge) is torn down first.
    std::unique_ptr<bridge_type> bridge_;
    solver_storage_type solver_;
};

}

#endif
