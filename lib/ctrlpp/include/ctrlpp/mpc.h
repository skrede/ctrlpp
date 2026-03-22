#ifndef HPP_GUARD_CTRLPP_MPC_H
#define HPP_GUARD_CTRLPP_MPC_H

#include "ctrlpp/mpc/diagnostics.h"
#include "ctrlpp/mpc/qp_formulation.h"
#include "ctrlpp/mpc/qp_solver.h"
#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/terminal_set.h"

#include "ctrlpp/dare.h"
#include "ctrlpp/state_space.h"
#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <cstddef>
#include <optional>
#include <span>
#include <utility>
#include <variant>
#include <vector>

namespace ctrlpp {

template<typename Scalar, std::size_t NX, std::size_t NU>
struct mpc_config {
    int horizon{1};
    Matrix<Scalar, NX, NX> Q{Matrix<Scalar, NX, NX>::Identity()};
    Matrix<Scalar, NU, NU> R{Matrix<Scalar, NU, NU>::Identity()};
    std::optional<Matrix<Scalar, NX, NX>> Qf{};
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

template<typename Scalar, std::size_t NX, std::size_t NU, qp_solver Solver>
class mpc {
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);

public:
    mpc(const discrete_state_space<Scalar, NX, NU, NX>& system,
        const mpc_config<Scalar, NX, NU>& config)
        : config_{config}
        , system_{system}
        , u_prev_{Vector<Scalar, NU>::Zero()}
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

        // Terminal set constraint rows
        int n_terminal = 0;
        if (config_.terminal_constraint_set.has_value()) {
            n_terminal = std::visit([](const auto& s) -> int {
                using T = std::decay_t<decltype(s)>;
                if constexpr (std::is_same_v<T, ellipsoidal_set<Scalar, NX>>) {
                    return 2 * nx; // inner box approximation: 2 halfplanes per axis
                } else {
                    return static_cast<int>(s.H.rows());
                }
            }, config_.terminal_constraint_set.value());
        }

        n_con_ = (N + 1) * nx
            + (has_state_bounds ? N * nx : 0)
            + (has_input_bounds ? N * nu : 0)
            + (has_rate_bounds ? N * nu : 0)
            + n_terminal;

        // Compute terminal cost Qf
        Matrix<Scalar, NX, NX> Qf_actual;
        if (config_.Qf.has_value()) {
            Qf_actual = config_.Qf.value();
        } else {
            auto dare_result = dare<Scalar, NX, NU>(system_.A, system_.B, config_.Q, config_.R);
            Qf_actual = dare_result.value_or(config_.Q);
        }
        Qf_actual_ = Qf_actual;

        // Build QP structure (cold path)
        auto P = detail::build_cost_matrix<Scalar, NX, NU>(
            N, config_.Q, config_.R, Qf_actual, use_soft,
            config_.soft_penalty, config_.soft_state_penalty);
        auto A = detail::build_constraint_matrix<Scalar, NX, NU>(
            N, system_.A, system_.B,
            has_state_bounds, use_soft, has_input_bounds, has_rate_bounds,
            config_.terminal_constraint_set);

        Vector<Scalar, NX> x0_dummy = Vector<Scalar, NX>::Zero();
        auto [l, u] = detail::build_bounds_vectors<Scalar, NX, NU>(
            N, x0_dummy,
            config_.x_min, config_.x_max,
            config_.u_min, config_.u_max,
            config_.du_max, use_soft,
            config_.terminal_constraint_set);

        auto q = detail::build_cost_vector<Scalar, NX, NU>(N, n_dec_, config_.Q, Qf_actual);

        qp_problem<Scalar> problem{
            .P = std::move(P),
            .q = std::move(q),
            .A = std::move(A),
            .l = std::move(l),
            .u = std::move(u)
        };

        solver_.setup(problem);

        // Pre-allocate update vectors for hot path
        update_.q.resize(n_dec_);
        update_.q.setZero();
        update_.l.resize(n_con_);
        update_.l.setZero();
        update_.u.resize(n_con_);
        update_.u.setZero();
        // warm_x and warm_y intentionally left empty (size 0) until first successful solve
    }

    [[nodiscard]] auto solve(const Vector<Scalar, NX>& x0) -> std::optional<Vector<Scalar, NU>>
    {
        // Origin regulation: q is all zeros
        update_.q.setZero();
        return solve_impl(x0);
    }

    [[nodiscard]] auto solve(const Vector<Scalar, NX>& x0,
                             const Vector<Scalar, NX>& x_ref) -> std::optional<Vector<Scalar, NU>>
    {
        int N = config_.horizon;

        // Update q vector for setpoint tracking
        Vector<Scalar, NX> Qxr = config_.Q * x_ref;
        for (int k = 0; k < N; ++k)
            update_.q.segment(k * nx, nx) = -Qxr;
        update_.q.segment(N * nx, nx) = -(Qf_actual_ * x_ref);

        // Zero the input portion and slack portion
        update_.q.segment(n_x_total_, n_dec_ - n_x_total_).setZero();

        return solve_impl(x0);
    }

    [[nodiscard]] auto solve(const Vector<Scalar, NX>& x0,
                             std::span<const Vector<Scalar, NX>> x_ref) -> std::optional<Vector<Scalar, NU>>
    {
        int N = config_.horizon;

        // Update q vector for trajectory tracking
        for (int k = 0; k < N; ++k)
            update_.q.segment(k * nx, nx) = -(config_.Q * x_ref[static_cast<std::size_t>(k)]);
        update_.q.segment(N * nx, nx) = -(Qf_actual_ * x_ref[static_cast<std::size_t>(N)]);

        // Zero the input portion and slack portion
        update_.q.segment(n_x_total_, n_dec_ - n_x_total_).setZero();

        return solve_impl(x0);
    }

    [[nodiscard]] auto trajectory() const
        -> std::pair<std::vector<Vector<Scalar, NX>>, std::vector<Vector<Scalar, NU>>>
    {
        int N = config_.horizon;
        std::vector<Vector<Scalar, NX>> states;
        std::vector<Vector<Scalar, NU>> inputs;
        states.reserve(static_cast<std::size_t>(N + 1));
        inputs.reserve(static_cast<std::size_t>(N));

        for (int k = 0; k <= N; ++k)
            states.push_back(last_primal_.segment(k * nx, nx));
        for (int k = 0; k < N; ++k)
            inputs.push_back(last_primal_.segment(n_x_total_ + k * nu, nu));

        return {std::move(states), std::move(inputs)};
    }

    [[nodiscard]] auto diagnostics() const -> mpc_diagnostics<Scalar>
    {
        return last_diagnostics_;
    }

private:
    [[nodiscard]] auto solve_impl(const Vector<Scalar, NX>& x0) -> std::optional<Vector<Scalar, NU>>
    {
        int N = config_.horizon;

        // Rebuild bounds for this solve
        auto [l, u] = detail::build_bounds_vectors<Scalar, NX, NU>(
            N, x0,
            config_.x_min, config_.x_max,
            config_.u_min, config_.u_max,
            config_.du_max,
            (config_.x_min.has_value() || config_.x_max.has_value())
                && !config_.hard_state_constraints,
            config_.terminal_constraint_set);
        update_.l = std::move(l);
        update_.u = std::move(u);

        // Update rate constraint bounds for k=0 using u_prev_
        if (config_.du_max.has_value()) {
            int n_dyn = (N + 1) * nx;
            bool has_state_bounds = config_.x_min.has_value() || config_.x_max.has_value();
            bool has_input_bounds = config_.u_min.has_value() || config_.u_max.has_value();
            int rate_row = n_dyn
                + (has_state_bounds ? N * nx : 0)
                + (has_input_bounds ? N * nu : 0);

            // k=0: constraint is u_0 - u_prev, so bounds shift by u_prev
            update_.l.segment(rate_row, nu) = -config_.du_max.value() + u_prev_;
            update_.u.segment(rate_row, nu) = config_.du_max.value() + u_prev_;
        }

        // Warm-start from previous solution
        if (has_solution_) {
            update_.warm_x = last_primal_;
            update_.warm_y = last_dual_;
        }

        auto result = solver_.solve(update_);

        // Populate diagnostics
        last_diagnostics_ = mpc_diagnostics<Scalar>{
            .status = result.status,
            .iterations = result.iterations,
            .solve_time = result.solve_time,
            .cost = result.objective,
            .primal_residual = result.primal_residual,
            .dual_residual = result.dual_residual,
            .max_constraint_violation = Scalar{0}
        };

        if (result.status != solve_status::optimal
            && result.status != solve_status::solved_inaccurate) {
            return std::nullopt;
        }

        // Store solution for warm-starting
        last_primal_ = std::move(result.x);
        last_dual_ = std::move(result.y);
        has_solution_ = true;

        // Extract u_0
        Vector<Scalar, NU> u0 = last_primal_.segment(n_x_total_, nu);
        u_prev_ = u0;

        return u0;
    }

    Solver solver_{};
    mpc_config<Scalar, NX, NU> config_;
    discrete_state_space<Scalar, NX, NU, NX> system_;
    Matrix<Scalar, NX, NX> Qf_actual_;

    qp_update<Scalar> update_{};
    Eigen::VectorX<Scalar> last_primal_{};
    Eigen::VectorX<Scalar> last_dual_{};
    mpc_diagnostics<Scalar> last_diagnostics_{};
    Vector<Scalar, NU> u_prev_;
    bool has_solution_{false};

    int n_dec_{};
    int n_con_{};
    int n_x_total_{};
    int n_u_total_{};
};

}

#endif
