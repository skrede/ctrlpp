#ifndef HPP_GUARD_CTRLPP_MPC_H
#define HPP_GUARD_CTRLPP_MPC_H

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/control/dare.h"

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/qp_solver.h"
#include "ctrlpp/mpc/diagnostics.h"
#include "ctrlpp/mpc/terminal_set.h"
#include "ctrlpp/mpc/qp_formulation.h"

#include "ctrlpp/model/state_space.h"

#include <Eigen/Dense>

#include <span>
#include <limits>
#include <vector>
#include <cstddef>
#include <utility>
#include <variant>
#include <optional>

namespace ctrlpp
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY = NX>
struct mpc_config
{
    int horizon{1};
    Matrix<Scalar, NY, NY> Q{Matrix<Scalar, NY, NY>::Identity()};
    Matrix<Scalar, NU, NU> R{Matrix<Scalar, NU, NU>::Identity()};
    std::optional<Matrix<Scalar, NY, NY>> Qf{};
    std::optional<Vector<Scalar, NU>> u_min{};
    std::optional<Vector<Scalar, NU>> u_max{};
    std::optional<Vector<Scalar, NX>> x_min{};
    std::optional<Vector<Scalar, NX>> x_max{};
    std::optional<Vector<Scalar, NU>> du_max{};
    Scalar soft_penalty{Scalar{1e4}};
    std::optional<Vector<Scalar, NX>> soft_state_penalty{};
    std::optional<terminal_set<Scalar, NX>> terminal_constraint_set{};
    bool hard_state_constraints{false};
};

template <typename Scalar, std::size_t NX, std::size_t NU, qp_solver Solver, std::size_t NY = NX>
class mpc
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);
    static constexpr int ny = static_cast<int>(NY);

public:
    /// @brief Validating factory; the only construction path on an
    /// exception-free build. Builds the solver with its own defaults and chains
    /// into the solver-taking overload below.
    ///
    /// Rejections, checked in order:
    ///  * horizon <= 0                   -> controller_construction_error::non_positive_horizon
    ///  * horizon above the representable
    ///    bound of the derived dimensions -> controller_construction_error::horizon_overflow
    static auto create(const discrete_state_space<Scalar, NX, NU, NY>& system, const mpc_config<Scalar, NX, NU, NY>& config)
        -> expected<mpc, controller_construction_error>
    {
        return create(system, config, Solver{});
    }

    /// @brief Validating factory taking a caller-supplied, pre-configured
    /// solver, e.g. `mpc<...>::create(sys, cfg, osqp_solver{qp_preset::speed})`
    /// to skip per-step polishing on the warm-resolve MPC path. The solver is
    /// moved in before the initial QP is posed, so its settings govern setup.
    ///
    /// The horizon is the only runtime quantity that scales the posed problem,
    /// and it is validated here, before any dimension product is formed and
    /// before any storage is reserved. Rejections, checked in order:
    ///  * horizon <= 0                     -> controller_construction_error::non_positive_horizon
    ///  * horizon > horizon_bound(config)   -> controller_construction_error::horizon_overflow
    ///
    /// The overflow bound is a representability condition on the horizon's own
    /// type, not a chosen ceiling. At the worst-case configuration (soft state
    /// bounds, input bounds and rate bounds all present) the horizon N scales
    /// two dimensions:
    ///   decision vector : (N+1)*nx + N*nu + N*nx  = N*(2*nx + nu) + nx
    ///   constraint rows : (N+1)*nx + N*nx + 2*N*nu
    ///                                             = N*(2*nx + 2*nu) + nx + n_terminal
    /// so the largest per-step contribution is 2*nx + 2*nu and the
    /// horizon-independent part is nx + n_terminal (the terminal rows come from
    /// the caller's terminal set and do not scale with the horizon). Both
    /// products therefore stay representable exactly when
    ///   horizon <= (max<int> - (nx + n_terminal)) / (2*nx + 2*nu).
    ///
    /// Two alternatives are deliberately not implemented. Validating at the
    /// first solve would surface a configuration error at the first control
    /// step, the worst possible moment. Clamping the horizon to one would turn a
    /// caller mistake into a silently different controller.
    static auto create(const discrete_state_space<Scalar, NX, NU, NY>& system, const mpc_config<Scalar, NX, NU, NY>& config, Solver solver)
        -> expected<mpc, controller_construction_error>
    {
        if(config.horizon <= 0)
            return unexpected(controller_construction_error::non_positive_horizon);
        if(config.horizon > horizon_bound(config))
            return unexpected(controller_construction_error::horizon_overflow);

        return mpc{unchecked_t{}, system, config, std::move(solver)};
    }

    // Unified soft-constraint / failure contract (shared by mpc and nmpc).
    //
    // solve() returns ctrlpp::expected<solve_output<Scalar, NU>, solver_error>.
    //   * SUCCESS branch: a usable control input plus a soft solve_result_status.
    //     - optimal            -> converged
    //     - solved_inaccurate  -> solved_inaccurate
    //     - max_iterations     -> budget_exhausted (best iterate still returned)
    //     - time_limit         -> budget_exhausted (best iterate still returned)
    //     The applied input is reached explicitly through `->input`; there is no
    //     implicit conversion to Vector, so the soft status can never be dropped.
    //   * ERROR branch: a hard failure with NO input.
    //     - infeasible                     -> infeasible
    //     - unbounded / non_convex / error -> invalid_problem
    //     - one-time solver setup failed   -> setup_incomplete
    //     - result too short for the problem -> invalid_backend_result
    //     On the error branch the controller does NOT update its internal u_prev
    //     and applies no hidden fallback input, so a failed solve never warms the
    //     rate constraints from a phantom input. Use set_applied_input to record
    //     the input the caller actually commanded.
    auto solve(const Vector<Scalar, NX>& x0) -> expected<solve_output<Scalar, NU>, solver_error>
    {
        update_.q.setZero();
        return solve_impl(x0);
    }

    auto solve(const Vector<Scalar, NX>& x0, const Vector<Scalar, NY>& y_ref) -> expected<solve_output<Scalar, NU>, solver_error>
    {
        int N = config_.horizon;
        Vector<Scalar, NX> CtQ_yref = CtQ_ * y_ref;
        for(int k = 0; k < N; ++k)
            update_.q.segment(k * nx, nx) = -CtQ_yref;
        update_.q.segment(N * nx, nx) = -(Qf_linear_ * y_ref);
        update_.q.segment(n_x_total_, n_dec_ - n_x_total_).setZero();
        return solve_impl(x0);
    }

    auto solve(const Vector<Scalar, NX>& x0, std::span<const Vector<Scalar, NY>> y_ref) -> expected<solve_output<Scalar, NU>, solver_error>
    {
        int N = config_.horizon;
        // The tracking overload reads y_ref[0..N] (N+1 references). An undersized
        // span would otherwise overrun; reject it via the error branch.
        if(y_ref.size() < static_cast<std::size_t>(N + 1))
            return unexpected(solver_error::invalid_problem);
        for(int k = 0; k < N; ++k)
            update_.q.segment(k * nx, nx) = -(CtQ_ * y_ref[static_cast<std::size_t>(k)]);
        update_.q.segment(N * nx, nx) = -(Qf_linear_ * y_ref[static_cast<std::size_t>(N)]);
        update_.q.segment(n_x_total_, n_dec_ - n_x_total_).setZero();
        return solve_impl(x0);
    }

    /// @brief Record the control input the caller actually commanded.
    ///
    /// The internal u_prev (which anchors the rate constraint) is updated only on
    /// a successful solve. When the caller overrides the commanded input (for
    /// example after an error branch, or when saturating externally), this lets it
    /// keep the rate-constraint reference consistent with what was truly applied.
    void set_applied_input(const Vector<Scalar, NU>& u) { u_prev_ = u; }

    // Guarded: returns the error branch before the first valid solve, so a caller
    // can never read stale or default-initialized primal data.
    auto trajectory() const -> expected<std::pair<std::vector<Vector<Scalar, NX>>, std::vector<Vector<Scalar, NU>>>, solver_error>
    {
        if(!has_solution_)
            return unexpected(solver_error::setup_incomplete);

        // The slices below need no length check of their own: has_solution_ is
        // set only by extract_solution, which stores the primal only after
        // checking that it covers the decision dimension these offsets are
        // derived from. A stored primal is therefore long enough by construction.
        int N = config_.horizon;
        std::vector<Vector<Scalar, NX>> states;
        std::vector<Vector<Scalar, NU>> inputs;
        states.reserve(static_cast<std::size_t>(N + 1));
        inputs.reserve(static_cast<std::size_t>(N));

        for(int k = 0; k <= N; ++k)
            states.push_back(last_primal_.segment(k * nx, nx));
        for(int k = 0; k < N; ++k)
            inputs.push_back(last_primal_.segment(n_x_total_ + k * nu, nu));

        return std::pair{std::move(states), std::move(inputs)};
    }

    auto diagnostics() const -> mpc_diagnostics<Scalar> { return last_diagnostics_; }

private:
    /// @brief Tag selecting the non-validating constructor reserved for `create`.
    struct unchecked_t
    {
        explicit unchecked_t() = default;
    };

    /// @brief Construct from a configuration already validated by `create`.
    mpc(unchecked_t, const discrete_state_space<Scalar, NX, NU, NY>& system, const mpc_config<Scalar, NX, NU, NY>& config, Solver solver) : solver_{std::move(solver)}, config_{config}, system_{system}, u_prev_{Vector<Scalar, NU>::Zero()}
    {
        precompute_output_weights();
        compute_dimensions();
        compute_terminal_cost();
        build_initial_qp();
        allocate_update_vectors();
    }

    /// @brief Largest horizon whose derived decision and constraint dimensions
    /// are still representable in the horizon's own type. See `create` for
    /// the derivation; this forms no product of its own.
    static auto horizon_bound(const mpc_config<Scalar, NX, NU, NY>& config) -> int
    {
        constexpr int per_step = 2 * nx + 2 * nu;
        const int constant_dimensions = nx + detail::terminal_constraint_rows<Scalar, NX>(config.terminal_constraint_set);
        return (std::numeric_limits<int>::max() - constant_dimensions) / per_step;
    }

    /// @brief Derive the decision width, the constraint height, and the row the
    /// rate block starts at from the configuration's own optionals.
    ///
    /// The same derivation sizes the constraint matrix and the bounds vectors, so
    /// taking all of it from one helper is what keeps the artifacts handed to the
    /// backend the same height as the problem this controller believes it posed.
    void compute_dimensions()
    {
        int N = config_.horizon;
        n_x_total_ = (N + 1) * nx;
        n_u_total_ = N * nu;

        n_dec_ = n_x_total_ + n_u_total_ + detail::slack_columns<NX>(N, uses_soft_state_constraints());

        auto const layout = constraint_layout();
        n_con_ = layout.n_con;
        rate_row_ = layout.rate_row;
    }

    /// @brief Whether the state bounds are softened with slack variables.
    ///
    /// Not derived from an optional: state bounds being present is a necessary
    /// condition, but enforcing them hard instead is the caller's policy choice.
    auto uses_soft_state_constraints() const -> bool
    {
        return (config_.x_min.has_value() || config_.x_max.has_value()) && !config_.hard_state_constraints;
    }

    auto constraint_layout() const -> detail::qp_constraint_layout
    {
        return detail::constraint_layout<Scalar, NX, NU>(config_.horizon, config_.x_min, config_.x_max, config_.u_min, config_.u_max, config_.du_max, config_.terminal_constraint_set);
    }

    void precompute_output_weights()
    {
        Q_state_ = system_.C.transpose() * config_.Q * system_.C;
        CtQ_ = system_.C.transpose() * config_.Q;
    }

    void compute_terminal_cost()
    {
        if(config_.Qf.has_value())
        {
            Qf_state_ = system_.C.transpose() * (*config_.Qf) * system_.C;
            Qf_linear_ = system_.C.transpose() * (*config_.Qf);
        }
        else
        {
            auto dare_result = dare<Scalar, NX, NU>(system_.A, system_.B, Q_state_, config_.R);
            Qf_state_ = dare_result ? dare_result->P : Q_state_;
            // Map state-space Qf back to output space for the linear tracking term:
            // linear term = -Qf_state * C_pinv * y_ref, where C_pinv = C' * (C*C')^{-1}
            Matrix<Scalar, NY, NY> CCt = system_.C * system_.C.transpose();
            Qf_linear_ = Qf_state_ * system_.C.transpose() * CCt.colPivHouseholderQr().solve(Matrix<Scalar, NY, NY>::Identity());
        }
    }

    void build_initial_qp()
    {
        int N = config_.horizon;
        bool use_soft = uses_soft_state_constraints();

        // Both builders read the same optionals, so the matrix and the bounds
        // vectors cannot come out of this function with different row counts.
        auto P = detail::build_cost_matrix<Scalar, NX, NU>(N, Q_state_, config_.R, Qf_state_, use_soft, config_.soft_penalty, config_.soft_state_penalty);
        auto A = detail::build_constraint_matrix<Scalar, NX, NU>(N, system_.A, system_.B, config_.x_min, config_.x_max, use_soft, config_.u_min, config_.u_max, config_.du_max, config_.terminal_constraint_set);

        Vector<Scalar, NX> x0_dummy = Vector<Scalar, NX>::Zero();
        auto [l, u] = detail::build_bounds_vectors<Scalar, NX, NU>(N, x0_dummy, config_.x_min, config_.x_max, config_.u_min, config_.u_max, config_.du_max, config_.terminal_constraint_set);

        auto q = detail::build_cost_vector<Scalar, NX, NU>(N, n_dec_, Q_state_, Qf_state_);

        qp_problem<Scalar> problem{.P = std::move(P), .q = std::move(q), .A = std::move(A), .l = std::move(l), .u = std::move(u)};
        setup_failed_ = !detail::setup_qp_solver(solver_, problem).has_value();
    }

    void allocate_update_vectors()
    {
        update_.q.resize(n_dec_);
        update_.q.setZero();
        update_.l.resize(n_con_);
        update_.l.setZero();
        update_.u.resize(n_con_);
        update_.u.setZero();
    }

    auto solve_impl(const Vector<Scalar, NX>& x0) -> expected<solve_output<Scalar, NU>, solver_error>
    {
        if(setup_failed_)
            return unexpected(solver_error::setup_incomplete);

        rebuild_bounds(x0);
        apply_rate_constraint_update();
        apply_warm_start();

        auto result = solver_.solve(update_);

        populate_diagnostics(result);

        // WIDENED accept-set: budget-limited iterates (max_iterations/time_limit)
        // now reach the caller on the SUCCESS branch tagged budget_exhausted.
        switch(result.status)
        {
        case solve_status::optimal:
            return extract_solution(result, solve_result_status::converged);
        case solve_status::solved_inaccurate:
            return extract_solution(result, solve_result_status::solved_inaccurate);
        case solve_status::max_iterations:
        case solve_status::time_limit:
            return extract_solution(result, solve_result_status::budget_exhausted);
        case solve_status::infeasible:
            return unexpected(solver_error::infeasible);
        case solve_status::unbounded:
        case solve_status::non_convex:
        case solve_status::error:
        default:
            return unexpected(solver_error::invalid_problem);
        }
    }

    void rebuild_bounds(const Vector<Scalar, NX>& x0)
    {
        int N = config_.horizon;
        auto [l, u] = detail::build_bounds_vectors<Scalar, NX, NU>(N,
                                                                   x0,
                                                                   config_.x_min,
                                                                   config_.x_max,
                                                                   config_.u_min,
                                                                   config_.u_max,
                                                                   config_.du_max,
                                                                   config_.terminal_constraint_set);
        update_.l = std::move(l);
        update_.u = std::move(u);
    }

    void apply_rate_constraint_update()
    {
        if(!config_.du_max.has_value())
            return;

        // The offset comes from the same derivation that sized the bounds vectors
        // this function writes into, so it cannot address a different block.
        update_.l.segment(rate_row_, nu) = -(*config_.du_max) + u_prev_;
        update_.u.segment(rate_row_, nu) = *config_.du_max + u_prev_;
    }

    void apply_warm_start()
    {
        if(has_solution_)
        {
            update_.warm_x = last_primal_;
            update_.warm_y = last_dual_;
        }
    }

    void populate_diagnostics(const qp_result<Scalar>& result)
    {
        last_diagnostics_ = mpc_diagnostics<Scalar>{.status = result.status,
                                                    .iterations = result.iterations,
                                                    .solve_time = result.solve_time,
                                                    .cost = result.objective,
                                                    .primal_residual = result.primal_residual,
                                                    .dual_residual = result.dual_residual,
                                                    .max_constraint_violation = Scalar{0}};
    }

    /// @brief Consume an accepted backend result: validate its reported shape,
    /// then extract the applied input.
    ///
    /// A status is not a shape. The accept-set above is decided purely from what
    /// the backend reports, and a backend that reports an accepted status can
    /// still return a primal shorter than the decision dimension or a dual
    /// shorter than the constraint count. The slice taken below sits at an offset
    /// derived from the horizon, so a short primal makes it read past the end of
    /// the backend's own storage, and a short dual is handed straight back as the
    /// next warm start. Both reported lengths are therefore compared against the
    /// dimensions computed at construction BEFORE either vector is moved from or
    /// indexed, and a violation leaves every member untouched.
    ///
    /// A longer-than-required result is not rejected: it is readable, and how
    /// much storage a backend returns beyond the posed problem is its own affair.
    /// The condition checked here is exactly the one that makes the reads legal.
    auto extract_solution(qp_result<Scalar>& result, solve_result_status status) -> expected<solve_output<Scalar, NU>, solver_error>
    {
        if(result.x.size() < static_cast<Eigen::Index>(n_dec_) || result.y.size() < static_cast<Eigen::Index>(n_con_))
        {
            // Correct the diagnostics the backend's own status just populated, so
            // a reader there is not told the solve was optimal when its answer
            // was discarded.
            last_diagnostics_.status = solve_status::invalid_backend_result;
            return unexpected(solver_error::invalid_backend_result);
        }

        last_primal_ = std::move(result.x);
        last_dual_ = std::move(result.y);
        has_solution_ = true;

        // The internal u_prev update lives on the SUCCESS branch only (this method
        // is reached solely from a successful solve), so a failed solve never
        // warms the rate constraint from a phantom input.
        Vector<Scalar, NU> u0 = last_primal_.segment(n_x_total_, nu);
        u_prev_ = u0;
        return solve_output<Scalar, NU>{.input = std::move(u0), .status = status};
    }

    Solver solver_{};
    mpc_config<Scalar, NX, NU, NY> config_;
    discrete_state_space<Scalar, NX, NU, NY> system_;
    Matrix<Scalar, NX, NX> Q_state_{};
    Matrix<Scalar, NX, NX> Qf_state_{};
    Matrix<Scalar, NX, NY> CtQ_{};
    Matrix<Scalar, NX, NY> Qf_linear_{};

    qp_update<Scalar> update_{};
    Eigen::VectorX<Scalar> last_primal_{};
    Eigen::VectorX<Scalar> last_dual_{};
    mpc_diagnostics<Scalar> last_diagnostics_{};
    Vector<Scalar, NU> u_prev_;
    bool has_solution_{false};
    bool setup_failed_{false};

    int n_dec_{};
    int n_con_{};
    int n_x_total_{};
    int n_u_total_{};
    int rate_row_{};
};

}

#endif
