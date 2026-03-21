#ifndef HPP_GUARD_CTRLPP_NMPC_H
#define HPP_GUARD_CTRLPP_NMPC_H

#include "ctrlpp/mpc/diagnostics.h"
#include "ctrlpp/mpc/dynamics_model.h"
#include "ctrlpp/mpc/nlp_formulation.h"
#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/nmpc_config.h"
#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace ctrlpp {

template<typename Scalar, std::size_t NX, std::size_t NU,
         nlp_solver Solver,
         dynamics_model<Scalar, NX, NU> Dynamics>
class nmpc {
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);

public:
    nmpc(Dynamics dynamics, const nmpc_config<Scalar, NX, NU>& config)
        : dynamics_{std::move(dynamics)}
        , config_{config}
        , N_{config.horizon}
        , n_vars_{(N_ + 1) * nx + N_ * nu}
        , state_{std::make_shared<nmpc_formulation_state<Scalar, NX, NU>>()}
    {
        // Initialize reference trajectory to zeros (regulation)
        state_->x_ref.resize(static_cast<std::size_t>(N_ + 1),
                             Vector<Scalar, NX>::Zero());

        // Build NLP problem with formulation state captured in lambdas
        problem_ = detail::build_nmpc_problem<Scalar, NX, NU>(
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
        for (int k = 0; k < N_; ++k) {
            inputs.push_back(last_solution_.segment((N_ + 1) * nx + k * nu, nu));
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
            .max_constraint_violation = result.primal_residual
        };

        if (result.status != solve_status::optimal
            && result.status != solve_status::solved_inaccurate) {
            return std::nullopt;
        }

        // Store full solution
        last_solution_ = result.x;

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

        // Extract u_0
        Vector<Scalar, NU> u0 = last_solution_.segment(u_offset, nu);
        u_prev_ = u0;

        return u0;
    }

    Dynamics dynamics_;
    nmpc_config<Scalar, NX, NU> config_;
    int N_;
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
