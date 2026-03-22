#ifndef HPP_GUARD_CTRLPP_NMPC_H
#define HPP_GUARD_CTRLPP_NMPC_H

#include "ctrlpp/mpc/diagnostics.h"
#include "ctrlpp/mpc/dynamics_model.h"
#include "ctrlpp/mpc/nlp_formulation.h"
#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/nmpc_config.h"
#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace ctrlpp {

template<typename Scalar, std::size_t NX, std::size_t NU,
         nlp_solver Solver,
         dynamics_model<Scalar, NX, NU> Dynamics,
         std::size_t NC = 0, std::size_t NTC = 0>
class nmpc {
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);
    static constexpr int nc = static_cast<int>(NC);
    static constexpr int ntc = static_cast<int>(NTC);

public:
    nmpc(Dynamics dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config)
        : dynamics_{std::move(dynamics)}
        , config_{config}
        , N_{config.horizon}
        , has_path_slack_{config.soft_constraints && NC > 0 && config.path_constraint.has_value()}
        , has_term_slack_{config.soft_constraints && NTC > 0 && config.terminal_constraint.has_value()}
        , n_path_slack_{has_path_slack_ ? N_ * nc : 0}
        , n_term_slack_{has_term_slack_ ? ntc : 0}
        , n_vars_{(N_ + 1) * nx + N_ * nu + n_path_slack_ + n_term_slack_}
        , state_{std::make_shared<nmpc_formulation_state<Scalar, NX, NU>>()}
    {
        // Initialize reference trajectory to zeros (regulation)
        state_->x_ref.resize(static_cast<std::size_t>(N_ + 1),
                             Vector<Scalar, NX>::Zero());

        // Build NLP problem with formulation state captured in lambdas
        problem_ = detail::build_nmpc_problem<Scalar, NX, NU, NC, NTC>(
            dynamics_, config_, state_);

        solver_.setup(problem_);

        // Pre-allocate warm-start vector (zeros for first solve)
        warm_z_ = Eigen::VectorX<Scalar>::Zero(n_vars_);
    }

    [[nodiscard]] auto solve(const Vector<Scalar, NX>& x0)
        -> std::optional<Vector<Scalar, NU>>
    {
        // Regulation: reference is all zeros
        for (auto& ref : state_->x_ref) {
            ref.setZero();
        }
        return solve_impl(x0);
    }

    [[nodiscard]] auto solve(const Vector<Scalar, NX>& x0,
                              const Vector<Scalar, NX>& x_ref)
        -> std::optional<Vector<Scalar, NU>>
    {
        // Setpoint tracking: constant reference
        for (auto& ref : state_->x_ref) {
            ref = x_ref;
        }
        return solve_impl(x0);
    }

    [[nodiscard]] auto solve(const Vector<Scalar, NX>& x0,
                              std::span<const Vector<Scalar, NX>> x_ref)
        -> std::optional<Vector<Scalar, NU>>
    {
        // Trajectory tracking: copy reference sequence
        const auto len = std::min(x_ref.size(), state_->x_ref.size());
        for (std::size_t k = 0; k < len; ++k) {
            state_->x_ref[k] = x_ref[k];
        }
        // If span is shorter, repeat last element
        if (!x_ref.empty()) {
            for (std::size_t k = len; k < state_->x_ref.size(); ++k) {
                state_->x_ref[k] = x_ref.back();
            }
        }
        return solve_impl(x0);
    }

    [[nodiscard]] auto trajectory() const
        -> std::pair<std::vector<Vector<Scalar, NX>>, std::vector<Vector<Scalar, NU>>>
    {
        std::vector<Vector<Scalar, NX>> states;
        std::vector<Vector<Scalar, NU>> inputs;
        states.reserve(static_cast<std::size_t>(N_ + 1));
        inputs.reserve(static_cast<std::size_t>(N_));

        for (int k = 0; k <= N_; ++k) {
            states.push_back(last_solution_.segment(k * nx, nx));
        }
        const int u_offset = (N_ + 1) * nx;
        for (int k = 0; k < N_; ++k) {
            inputs.push_back(last_solution_.segment(u_offset + k * nu, nu));
        }

        return {std::move(states), std::move(inputs)};
    }

    [[nodiscard]] auto diagnostics() const -> mpc_diagnostics<Scalar>
    {
        return last_diagnostics_;
    }

private:
    [[nodiscard]] auto solve_impl(const Vector<Scalar, NX>& x0)
        -> std::optional<Vector<Scalar, NU>>
    {
        // Update formulation state (lambdas read from this shared state)
        state_->x0 = x0;
        state_->u_prev = u_prev_;

        // Build update with warm-start initial guess
        nlp_update<Scalar> update;
        update.x0 = warm_z_;

        auto result = solver_.solve(update);

        // Populate diagnostics
        last_diagnostics_ = mpc_diagnostics<Scalar>{
            .status = result.status,
            .iterations = result.iterations,
            .solve_time = result.solve_time,
            .cost = result.objective,
            .primal_residual = result.primal_residual,
            .dual_residual = Scalar{0},
            .max_constraint_violation = result.primal_residual,
            .max_path_constraint_violation = Scalar{0},
            .max_terminal_constraint_violation = Scalar{0},
            .total_slack = Scalar{0}
        };

        if (result.status != solve_status::optimal
            && result.status != solve_status::solved_inaccurate) {
            return std::nullopt;
        }

        // Store full solution
        last_solution_ = result.x;

        // Populate constraint diagnostics from solution
        populate_constraint_diagnostics(result.x);

        // Warm-start shift: move trajectory forward by one step
        warm_z_ = result.x;

        // Shift states: x[k] <- x[k+1] for k=0..N-1, duplicate last
        for (int k = 0; k < N_; ++k) {
            warm_z_.segment(k * nx, nx) = result.x.segment((k + 1) * nx, nx);
        }
        // Last state stays (already in place from copy)

        // Shift inputs: u[k] <- u[k+1] for k=0..N-2, duplicate last
        const int u_offset = (N_ + 1) * nx;
        for (int k = 0; k < N_ - 1; ++k) {
            warm_z_.segment(u_offset + k * nu, nu) =
                result.x.segment(u_offset + (k + 1) * nu, nu);
        }
        // Last input stays (already in place from copy)

        // Shift path slack: s[k] <- s[k+1] for k=0..N-2, zero last
        if (has_path_slack_) {
            const int slack_off = u_offset + N_ * nu;
            for (int k = 0; k < N_ - 1; ++k) {
                warm_z_.segment(slack_off + k * nc, nc) =
                    result.x.segment(slack_off + (k + 1) * nc, nc);
            }
            warm_z_.segment(slack_off + (N_ - 1) * nc, nc).setZero();
        }

        // Zero terminal slack
        if (has_term_slack_) {
            const int term_off = u_offset + N_ * nu + n_path_slack_;
            warm_z_.segment(term_off, ntc).setZero();
        }

        // Extract u_0
        Vector<Scalar, NU> u0 = last_solution_.segment(u_offset, nu);
        u_prev_ = u0;

        return u0;
    }

    void populate_constraint_diagnostics(const Eigen::VectorX<Scalar>& z)
    {
        const int u_offset = (N_ + 1) * nx;

        // Path constraint violations
        if constexpr (NC > 0) {
            if (config_.path_constraint) {
                const auto& g = *config_.path_constraint;
                Scalar max_viol{0};

                for (int k = 0; k < N_; ++k) {
                    Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + k * nx);
                    Eigen::Map<const Vector<Scalar, NU>> uk(z.data() + u_offset + k * nu);
                    Vector<Scalar, NC> gk = g(xk, uk);

                    for (int i = 0; i < nc; ++i) {
                        max_viol = std::max(max_viol, gk[i]);
                    }
                }

                last_diagnostics_.max_path_constraint_violation = max_viol;
            }
        }

        // Terminal constraint violations
        if constexpr (NTC > 0) {
            if (config_.terminal_constraint) {
                const auto& h = *config_.terminal_constraint;
                Eigen::Map<const Vector<Scalar, NX>> xN(z.data() + N_ * nx);
                Vector<Scalar, NTC> hN = h(xN);

                Scalar max_viol{0};
                for (int i = 0; i < ntc; ++i) {
                    max_viol = std::max(max_viol, hN[i]);
                }

                last_diagnostics_.max_terminal_constraint_violation = max_viol;
            }
        }

        // Total slack
        Scalar total{0};
        if (has_path_slack_) {
            const int slack_off = u_offset + N_ * nu;
            for (int k = 0; k < N_; ++k) {
                Eigen::Map<const Vector<Scalar, NC>> sk(z.data() + slack_off + k * nc);
                total += sk.sum();
            }
        }
        if (has_term_slack_) {
            const int term_off = u_offset + N_ * nu + n_path_slack_;
            Eigen::Map<const Vector<Scalar, NTC>> st(z.data() + term_off);
            total += st.sum();
        }
        last_diagnostics_.total_slack = total;
    }

    Dynamics dynamics_;
    nmpc_config<Scalar, NX, NU, NC, NTC> config_;
    int N_;
    bool has_path_slack_;
    bool has_term_slack_;
    int n_path_slack_;
    int n_term_slack_;
    int n_vars_;

    std::shared_ptr<nmpc_formulation_state<Scalar, NX, NU>> state_;
    nlp_problem<Scalar> problem_;
    Solver solver_{};

    Eigen::VectorX<Scalar> warm_z_;
    Eigen::VectorX<Scalar> last_solution_;
    mpc_diagnostics<Scalar> last_diagnostics_{};
    Vector<Scalar, NU> u_prev_{Vector<Scalar, NU>::Zero()};
};

}

#endif
