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
    mpc(const discrete_state_space<Scalar, NX, NU, NY>& system, const mpc_config<Scalar, NX, NU, NY>& config) : config_{config}, system_{system}, u_prev_{Vector<Scalar, NU>::Zero()}
    {
        precompute_output_weights();
        compute_dimensions();
        compute_terminal_cost();
        build_initial_qp();
        allocate_update_vectors();
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
    //     On the error branch the controller does NOT update its internal u_prev
    //     and applies no hidden fallback input, so a failed solve never warms the
    //     rate constraints from a phantom input. Use set_applied_input to record
    //     the input the caller actually commanded.
    [[nodiscard]] auto solve(const Vector<Scalar, NX>& x0) -> expected<solve_output<Scalar, NU>, solver_error>
    {
        update_.q.setZero();
        return solve_impl(x0);
    }

    [[nodiscard]] auto solve(const Vector<Scalar, NX>& x0, const Vector<Scalar, NY>& y_ref) -> expected<solve_output<Scalar, NU>, solver_error>
    {
        int N = config_.horizon;
        Vector<Scalar, NX> CtQ_yref = CtQ_ * y_ref;
        for(int k = 0; k < N; ++k)
            update_.q.segment(k * nx, nx) = -CtQ_yref;
        update_.q.segment(N * nx, nx) = -(Qf_linear_ * y_ref);
        update_.q.segment(n_x_total_, n_dec_ - n_x_total_).setZero();
        return solve_impl(x0);
    }

    [[nodiscard]] auto solve(const Vector<Scalar, NX>& x0, std::span<const Vector<Scalar, NY>> y_ref) -> expected<solve_output<Scalar, NU>, solver_error>
    {
        int N = config_.horizon;
        // The tracking overload reads y_ref[0..N] (N+1 references). An undersized
        // span would otherwise overrun; reject it via the error branch.
        if(y_ref.size() < static_cast<std::size_t>(N + 1))
            return unexpected<solver_error>{solver_error::invalid_problem};
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
    [[nodiscard]] auto trajectory() const -> expected<std::pair<std::vector<Vector<Scalar, NX>>, std::vector<Vector<Scalar, NU>>>, solver_error>
    {
        if(!has_solution_)
            return unexpected<solver_error>{solver_error::setup_incomplete};

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

    [[nodiscard]] auto diagnostics() const -> mpc_diagnostics<Scalar> { return last_diagnostics_; }

private:
    void compute_dimensions()
    {
        int N = config_.horizon;
        n_x_total_ = (N + 1) * nx;
        n_u_total_ = N * nu;

        bool has_state_bounds = config_.x_min.has_value() || config_.x_max.has_value();
        bool use_soft = has_state_bounds && !config_.hard_state_constraints;
        int n_slack = use_soft ? N * nx : 0;
        n_dec_ = n_x_total_ + n_u_total_ + n_slack;

        bool has_input_bounds = config_.u_min.has_value() || config_.u_max.has_value();
        bool has_rate_bounds = config_.du_max.has_value();
        int n_terminal = compute_terminal_constraint_rows();

        n_con_ = (N + 1) * nx + (has_state_bounds ? N * nx : 0) + (has_input_bounds ? N * nu : 0) + (has_rate_bounds ? N * nu : 0) + n_terminal;
    }

    int compute_terminal_constraint_rows() const
    {
        if(!config_.terminal_constraint_set.has_value())
            return 0;

        return std::visit(
            [](const auto& s) -> int
            {
                using T = std::decay_t<decltype(s)>;
                if constexpr(std::is_same_v<T, ellipsoidal_set<Scalar, NX>>)
                    return nx;
                else
                    return static_cast<int>(s.H.rows());
            },
            config_.terminal_constraint_set.value());
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
            Qf_state_ = system_.C.transpose() * config_.Qf.value() * system_.C;
            Qf_linear_ = system_.C.transpose() * config_.Qf.value();
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
        bool has_state_bounds = config_.x_min.has_value() || config_.x_max.has_value();
        bool use_soft = has_state_bounds && !config_.hard_state_constraints;

        auto P = detail::build_cost_matrix<Scalar, NX, NU>(N, Q_state_, config_.R, Qf_state_, use_soft, config_.soft_penalty, config_.soft_state_penalty);
        auto A = detail::build_constraint_matrix<Scalar, NX, NU>(N, system_.A, system_.B, has_state_bounds, use_soft, config_.u_min.has_value() || config_.u_max.has_value(), config_.du_max.has_value(), config_.terminal_constraint_set);

        Vector<Scalar, NX> x0_dummy = Vector<Scalar, NX>::Zero();
        auto [l, u] = detail::build_bounds_vectors<Scalar, NX, NU>(N, x0_dummy, config_.x_min, config_.x_max, config_.u_min, config_.u_max, config_.du_max, use_soft, config_.terminal_constraint_set);

        auto q = detail::build_cost_vector<Scalar, NX, NU>(N, n_dec_, Q_state_, Qf_state_);

        qp_problem<Scalar> problem{.P = std::move(P), .q = std::move(q), .A = std::move(A), .l = std::move(l), .u = std::move(u)};
        setup_failed_ = !detail::setup_qp_solver(solver_, problem);
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

    [[nodiscard]] auto solve_impl(const Vector<Scalar, NX>& x0) -> expected<solve_output<Scalar, NU>, solver_error>
    {
        if(setup_failed_)
            return unexpected<solver_error>{solver_error::setup_incomplete};

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
            return unexpected<solver_error>{solver_error::infeasible};
        case solve_status::unbounded:
        case solve_status::non_convex:
        case solve_status::error:
        default:
            return unexpected<solver_error>{solver_error::invalid_problem};
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
                                                                   (config_.x_min.has_value() || config_.x_max.has_value()) && !config_.hard_state_constraints,
                                                                   config_.terminal_constraint_set);
        update_.l = std::move(l);
        update_.u = std::move(u);
    }

    void apply_rate_constraint_update()
    {
        if(!config_.du_max.has_value())
            return;

        int N = config_.horizon;
        int n_dyn = (N + 1) * nx;
        bool has_state_bounds = config_.x_min.has_value() || config_.x_max.has_value();
        bool has_input_bounds = config_.u_min.has_value() || config_.u_max.has_value();
        int rate_row = n_dyn + (has_state_bounds ? N * nx : 0) + (has_input_bounds ? N * nu : 0);

        update_.l.segment(rate_row, nu) = -config_.du_max.value() + u_prev_;
        update_.u.segment(rate_row, nu) = config_.du_max.value() + u_prev_;
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

    auto extract_solution(qp_result<Scalar>& result, solve_result_status status) -> solve_output<Scalar, NU>
    {
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
};

}

#endif
