#ifndef HPP_GUARD_CTRLPP_MPC_NLOPT_SOLVER_H
#define HPP_GUARD_CTRLPP_MPC_NLOPT_SOLVER_H

/// @brief NLopt-backed NLP solver wrapper for nonlinear MPC and NMHE.
///
/// Setup fallibility is reported through `setup` returning
/// `ctrlpp::expected<void, nlopt_setup_error>`. Note that NLopt's C++ API throws
/// by upstream design (`setup` and `solve` contain every exception, standard or
/// not, and report it as an `nlopt_setup_error` or as `solve_status::error`,
/// except a forced unwind on libstdc++, which is re-raised), so this opt-in
/// backend requires exception support and is exempt from the library's embedded
/// no-exceptions floor.

#include "ctrlpp/expected.h"

#include "ctrlpp/mpc/nlp_types.h"
#include "ctrlpp/mpc/nlp_solver.h"

#include "ctrlpp/mpc/detail/nlopt_constraint_rows.h"

#include <nlopt.hpp>

#include <new>
#include <span>
#include <chrono>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <exception>
#include <type_traits>

#if defined(__GLIBCXX__)
#include <cxxabi.h>
#endif

namespace ctrlpp
{

/// NLopt algorithm selection for the NMPC solver backend.
/// Note: raw `mma` and raw `ccsaq` do not support equality constraints.
/// The `auglag_mma` and `auglag_ccsaq` variants wrap the inner solver
/// in NLOPT_AUGLAG_EQ so equality constraints are absorbed into the
/// outer penalty and inequalities pass through to the inner.
enum class nlopt_algorithm : std::uint8_t
{
    slsqp,
    mma,
    cobyla,
    isres,
    auglag_mma,
    auglag_ccsaq,
    ccsaq
};

template <typename Scalar>
struct nlopt_settings
{
    nlopt_algorithm algorithm{nlopt_algorithm::slsqp};
    Scalar ftol_rel{Scalar{1e-6}};
    Scalar xtol_rel{Scalar{1e-6}};
    int max_eval{500};
    Scalar max_time{Scalar{0}};
    Scalar constraint_tol{Scalar{1e-8}};
};

template <typename Scalar>
class nlopt_solver
{
    static_assert(std::is_same_v<Scalar, double>, "NLopt operates in double precision only");

public:
    using scalar_type = Scalar;

    explicit nlopt_solver(nlopt_settings<Scalar> settings = {}) : rows_{}, problem_{nullptr}, settings_{settings}, opt_{}, ready_{false} {}

    // Every NLopt callback is registered with the ADDRESS of this solver or of its constraint rows, so a
    // bitwise copy or move would leave the new optimizer calling back into the source, which every fallible
    // factory's return by value would reach. The copy and the move therefore re-run setup against the same
    // referenced problem, which the caller keeps at a stable address for the lifetime of the solver.
    nlopt_solver(const nlopt_solver& other) : rows_{}, problem_{other.problem_}, settings_{other.settings_}, opt_{}, ready_{false} { rebind_callbacks(); }

    nlopt_solver(nlopt_solver&& other) : rows_{}, problem_{other.problem_}, settings_{other.settings_}, opt_{}, ready_{false} { rebind_callbacks(); }

    auto operator=(const nlopt_solver& other) -> nlopt_solver&
    {
        if(this != &other)
        {
            settings_ = other.settings_;
            problem_ = other.problem_;
            rebind_callbacks();
        }
        return *this;
    }

    auto operator=(nlopt_solver&& other) -> nlopt_solver&
    {
        if(this != &other)
        {
            settings_ = other.settings_;
            problem_ = other.problem_;
            rebind_callbacks();
        }
        return *this;
    }

    ~nlopt_solver() = default;

    /// @brief Fallible setup: configures the NLopt optimizer from the problem
    /// definition. Returns an empty expected on success,
    /// `nlopt_setup_error::incompatible_equality_constraints` when the selected
    /// algorithm rejects the constraint structure,
    /// `nlopt_setup_error::allocation_failed` when an allocation fails while the
    /// optimizer is built or configured, and `nlopt_setup_error::backend_rejected`
    /// for any other error raised meanwhile, standard or not. A failed setup
    /// leaves the solver unusable: `solve` reports `solve_status::error` until a
    /// setup succeeds. This is the only setup shape.
    auto setup(const nlp_problem<Scalar>& problem) -> ctrlpp::expected<void, nlopt_setup_error>
    {
        using outcome = ctrlpp::expected<void, nlopt_setup_error>;
        ready_ = false;
        problem_ = &problem;
        const outcome configured = contained<outcome>([&] { return configure(problem); }, ctrlpp::unexpected(nlopt_setup_error::allocation_failed),
                                                      ctrlpp::unexpected(nlopt_setup_error::backend_rejected));
        ready_ = configured.has_value();
        return configured;
    }

    auto solve(const nlp_update<Scalar>& update) -> nlp_result<Scalar>
    {
        const auto t0 = std::chrono::steady_clock::now();
        nlp_result<Scalar> result = ready_ ? guarded_solve(update) : failed_result();
        result.solve_time = static_cast<Scalar>(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        return result;
    }

private:
    detail::nlopt_constraint_rows<Scalar> rows_;
    const nlp_problem<Scalar>* problem_;
    nlopt_settings<Scalar> settings_;
    nlopt::opt opt_;
    bool ready_;

    static auto to_nlopt_algorithm(nlopt_algorithm alg) -> nlopt::algorithm
    {
        switch(alg)
        {
        case nlopt_algorithm::slsqp:
            return nlopt::LD_SLSQP;
        case nlopt_algorithm::mma:
            return nlopt::LD_MMA;
        case nlopt_algorithm::cobyla:
            return nlopt::LN_COBYLA;
        case nlopt_algorithm::isres:
            return nlopt::GN_ISRES;
        case nlopt_algorithm::auglag_mma:
            return nlopt::AUGLAG_EQ;
        case nlopt_algorithm::auglag_ccsaq:
            return nlopt::AUGLAG_EQ;
        case nlopt_algorithm::ccsaq:
            return nlopt::LD_CCSAQ;
        }
        return nlopt::LD_SLSQP;
    }

    static auto to_nlopt_inner(nlopt_algorithm alg) -> nlopt::algorithm
    {
        switch(alg)
        {
        case nlopt_algorithm::auglag_mma:
            return nlopt::LD_MMA;
        case nlopt_algorithm::auglag_ccsaq:
            return nlopt::LD_CCSAQ;
        default:
            return nlopt::LD_SLSQP;
        }
    }

    auto configure(const nlp_problem<Scalar>& problem) -> ctrlpp::expected<void, nlopt_setup_error>
    {
        const bool rejects_equalities = settings_.algorithm == nlopt_algorithm::mma || settings_.algorithm == nlopt_algorithm::ccsaq;
        if(auto partitioned = rows_.partition(problem, rejects_equalities); !partitioned.has_value())
            return partitioned;

        opt_ = nlopt::opt(to_nlopt_algorithm(settings_.algorithm), static_cast<unsigned>(problem.n_vars));
        configure_variable_bounds(problem);
        opt_.set_min_objective(objective_callback, this);
        rows_.register_callbacks(opt_, static_cast<double>(settings_.constraint_tol));
        configure_stopping_criteria();
        if(settings_.algorithm == nlopt_algorithm::auglag_mma || settings_.algorithm == nlopt_algorithm::auglag_ccsaq)
            configure_local_optimizer(problem);
        return {};
    }

    void configure_local_optimizer(const nlp_problem<Scalar>& problem)
    {
        nlopt::opt local(to_nlopt_inner(settings_.algorithm), static_cast<unsigned>(problem.n_vars));
        local.set_ftol_rel(static_cast<double>(settings_.ftol_rel));
        local.set_xtol_rel(static_cast<double>(settings_.xtol_rel));
        local.set_maxeval(settings_.max_eval);
        if(settings_.max_time > Scalar{0})
            local.set_maxtime(static_cast<double>(settings_.max_time));
        opt_.set_local_optimizer(local);
    }

    // A forced unwind (glibc thread cancellation) is re-raised; one raised inside optimize is still ended by NLopt's callback wrapper.
    // https://gcc.gnu.org/onlinedocs/libstdc++/manual/using_exceptions.html#intro.using.exception.compat
    template <typename Result, typename Fn>
    static auto contained(Fn&& fn, Result on_allocation_failure, Result on_other_failure) -> Result
    {
        try
        {
            return std::forward<Fn>(fn)();
        }
#if defined(__GLIBCXX__)
        catch(const abi::__forced_unwind&)
        {
            throw;
        }
#endif
        catch(const std::bad_alloc&)
        {
            return on_allocation_failure;
        }
        catch(const std::exception&)
        {
            return on_other_failure;
        }
        catch(...)
        {
            return on_other_failure;
        }
    }

    // solve reports every failure through nlp_result::status, so no exception may cross it; run_optimizer's
    // round-off cost re-evaluation and the violation evaluation in solved_result run inside this guard.
    auto guarded_solve(const nlp_update<Scalar>& update) -> nlp_result<Scalar>
    {
        return contained<nlp_result<Scalar>>(
            [&]
            {
                auto x = to_stdvec(update.x0);
                project_onto_bounds(x);
                double opt_f{};
                const solve_status status = run_optimizer(x, opt_f);
                return solved_result(x, opt_f, status);
            },
            failed_result(), failed_result());
    }

    // NLopt refuses a start outside the bounds, so a valid box that excludes the previous
    // iterate would otherwise be refused. The comparisons leave inverted or NaN bounds
    // unchanged for NLopt to refuse, which std::clamp's lower <= upper precondition would not.
    void project_onto_bounds(std::vector<double>& x) const
    {
        const auto n = static_cast<Eigen::Index>(x.size());
        const bool has_lower = problem_->x_lower.size() == n;
        const bool has_upper = problem_->x_upper.size() == n;
        for(Eigen::Index i = 0; i < n; ++i)
        {
            double& xi = x[static_cast<std::size_t>(i)];
            if(has_upper && xi > static_cast<double>(problem_->x_upper[i]))
                xi = static_cast<double>(problem_->x_upper[i]);
            if(has_lower && xi < static_cast<double>(problem_->x_lower[i]))
                xi = static_cast<double>(problem_->x_lower[i]);
        }
    }

    auto run_optimizer(std::vector<double>& x, double& opt_f) -> solve_status
    {
        try
        {
            return map_result(opt_.optimize(x, opt_f));
        }
        catch(const nlopt::roundoff_limited&)
        {
            opt_f = problem_->cost_cb({x.data(), x.size()});
            return solve_status::solved_inaccurate;
        }
    }

    auto solved_result(const std::vector<double>& x, double opt_f, solve_status status) const -> nlp_result<Scalar>
    {
        nlp_result<Scalar> result{};
        result.status = status;
        result.x.resize(static_cast<Eigen::Index>(x.size()));
        for(Eigen::Index i = 0; i < result.x.size(); ++i)
            result.x[i] = static_cast<Scalar>(x[static_cast<std::size_t>(i)]);
        result.objective = static_cast<Scalar>(opt_f);
        result.iterations = opt_.get_numevals();
        result.primal_residual = rows_.max_violation(x);
        return result;
    }

    static auto failed_result() -> nlp_result<Scalar>
    {
        nlp_result<Scalar> result{};
        result.status = solve_status::error;
        return result;
    }

    void configure_variable_bounds(const nlp_problem<Scalar>& problem)
    {
        if(problem.x_lower.size() > 0)
            opt_.set_lower_bounds(to_stdvec(problem.x_lower));
        if(problem.x_upper.size() > 0)
            opt_.set_upper_bounds(to_stdvec(problem.x_upper));
    }

    // Re-register every callback against this object's address. The re-setup can
    // fail when the referenced problem changed since the source's setup; the copy
    // or move cannot report it, so it stays unready and its next solve reports
    // solve_status::error, as does a solver whose source was never set up.
    void rebind_callbacks()
    {
        ready_ = false;
        if(problem_ != nullptr)
            static_cast<void>(setup(*problem_));
    }

    void configure_stopping_criteria()
    {
        opt_.set_ftol_rel(static_cast<double>(settings_.ftol_rel));
        opt_.set_xtol_rel(static_cast<double>(settings_.xtol_rel));
        opt_.set_maxeval(settings_.max_eval);
        if(settings_.max_time > Scalar{0})
            opt_.set_maxtime(static_cast<double>(settings_.max_time));
    }

    static auto to_stdvec(const Eigen::VectorX<Scalar>& v) -> std::vector<double>
    {
        std::vector<double> out(static_cast<std::size_t>(v.size()));
        for(Eigen::Index i = 0; i < v.size(); ++i)
        {
            out[static_cast<std::size_t>(i)] = static_cast<double>(v[i]);
        }
        return out;
    }

    static constexpr auto map_result(nlopt::result rc) -> solve_status
    {
        switch(rc)
        {
        case nlopt::SUCCESS:
            [[fallthrough]];
        case nlopt::STOPVAL_REACHED:
            [[fallthrough]];
        case nlopt::FTOL_REACHED:
            [[fallthrough]];
        case nlopt::XTOL_REACHED:
            return solve_status::optimal;
        case nlopt::MAXEVAL_REACHED:
            return solve_status::max_iterations;
        case nlopt::MAXTIME_REACHED:
            return solve_status::time_limit;
        case nlopt::ROUNDOFF_LIMITED:
            return solve_status::solved_inaccurate;
        default:
            return solve_status::error;
        }
    }

    static auto objective_callback(const std::vector<double>& x, std::vector<double>& grad, void* data) -> double
    {
        auto* self = static_cast<nlopt_solver*>(data);
        std::span<const double> x_span{x.data(), x.size()};

        if(!grad.empty() && self->problem_->gradient_cb)
        {
            std::span<double> grad_span{grad.data(), grad.size()};
            self->problem_->gradient_cb(x_span, grad_span);
        }

        return self->problem_->cost_cb(x_span);
    }
};

}

#endif
