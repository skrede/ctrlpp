#ifndef HPP_GUARD_CTRLPP_MHE_MHE_NLP_FORMULATION_H
#define HPP_GUARD_CTRLPP_MHE_MHE_NLP_FORMULATION_H

/// @brief NLP formulation for nonlinear Moving Horizon Estimation with multiple shooting.
///
/// @cite diehl2009 -- Diehl et al., "Efficient Numerical Methods for Nonlinear MPC and Moving Horizon Estimation", 2009

#include "ctrlpp/mhe/mhe_config.h"

#include "ctrlpp/detail/numerical_diff.h"
#include "ctrlpp/model/differentiable_dynamics.h"
#include "ctrlpp/model/differentiable_measurement.h"
#include "ctrlpp/model/dynamics_model.h"
#include "ctrlpp/model/measurement_model.h"
#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <span>

namespace ctrlpp
{

/// Shared mutable state for the NMHE NLP formulation lambdas.
/// Updated by nmhe before each solve; lambdas capture a shared_ptr to this.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N>
struct nmhe_formulation_state
{
    Vector<Scalar, NX> arrival_state{Vector<Scalar, NX>::Zero()};
    Matrix<Scalar, NX, NX> arrival_P_inv{Matrix<Scalar, NX, NX>::Identity()};
    std::array<Vector<Scalar, NU>, N> u_window{};
    std::array<Vector<Scalar, NY>, N + 1> z_window{};
};

namespace detail
{

/// Build the NLP problem for nonlinear MHE with multiple shooting.
///
/// Decision vector layout:
///   z = [x_0, x_1, ..., x_N, s_path_0, ..., s_path_N, s_box_0, ..., s_box_N]
///   where path slack has (N+1)*NC entries (if NC > 0 and soft),
///   and box slack has (N+1)*NX entries (if box bounds present and soft).
///
/// Cost:
///   J = 0.5*(x_0 - x_bar)^T * (w * P_inv) * (x_0 - x_bar)
///     + sum_{k=0}^{N-1} 0.5 * w_k^T Q_inv w_k   where w_k = x_{k+1} - f(x_k, u_k)
///     + sum_{k=0}^{N}   0.5 * v_k^T R_inv v_k    where v_k = z_k - h(x_k)
///     + L1 slack penalties
///
/// Constraints:
///   Continuity:  x_{k+1} - f(x_k, u_k) = 0    for k = 0..N-1
///   Path:        g(x_k) - s_k <= 0              for k = 0..N (if NC > 0)
///   Residual:    |z_k,i - h_i(x_k)| <= bound_i  for k = 0..N (if residual_bound set)
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N, std::size_t NC = 0, typename Dynamics, typename Measurement>
    requires dynamics_model<Dynamics, Scalar, NX, NU> && measurement_model<Measurement, Scalar, NX, NY>
auto build_nmhe_problem(const Dynamics& dynamics,
                        const Measurement& measurement,
                        const nmhe_config<Scalar, NX, NU, NY, N, NC>& config,
                        std::shared_ptr<nmhe_formulation_state<Scalar, NX, NU, NY, N>> state,
                        const Matrix<Scalar, NX, NX>& Q_inv,
                        const Matrix<Scalar, NY, NY>& R_inv) -> nlp_problem<Scalar>
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int ny = static_cast<int>(NY);
    constexpr int nc = static_cast<int>(NC);
    constexpr int Ni = static_cast<int>(N);

    // Decision variable layout
    const bool has_path_constraint = NC > 0 && config.path_constraint.has_value();
    const bool has_path_slack = config.soft_constraints && has_path_constraint;
    const bool has_residual = config.residual_bound.has_value();

    const int n_states = (Ni + 1) * nx;
    const int n_path_slack = has_path_slack ? (Ni + 1) * nc : 0;
    const int n_vars = n_states + n_path_slack;

    const int x_offset = 0;
    const int path_slack_offset = n_states;

    // Constraint count
    const int n_continuity = Ni * nx;
    const int n_path_con = has_path_constraint ? (Ni + 1) * nc : 0;
    const int n_residual_con = has_residual ? (Ni + 1) * ny * 2 : 0;
    const int n_constraints = n_continuity + n_path_con + n_residual_con;

    const int continuity_start = 0;
    const int path_con_start = n_continuity;
    const int residual_con_start = path_con_start + n_path_con;

    const Scalar arrival_weight = config.arrival_cost_weight;
    // Cost function
    std::function<Scalar(std::span<const Scalar>)> cost = [=](std::span<const Scalar> z) -> Scalar
    {
        Scalar total{0};

        // Arrival cost: 0.5 * (x_0 - x_bar)^T * (w * P_inv) * (x_0 - x_bar)
        Eigen::Map<const Vector<Scalar, NX>> x0(z.data() + x_offset);
        Vector<Scalar, NX> dx0 = x0 - state->arrival_state;
        total += Scalar{0.5} * dx0.dot((arrival_weight * state->arrival_P_inv) * dx0);

        // Process noise cost: sum 0.5 * w_k^T Q_inv w_k
        for(int k = 0; k < Ni; ++k)
        {
            Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + x_offset + k * nx);
            Eigen::Map<const Vector<Scalar, NX>> xk1(z.data() + x_offset + (k + 1) * nx);

            Vector<Scalar, NX> wk = xk1 - dynamics(xk, state->u_window[static_cast<std::size_t>(k)]);
            total += Scalar{0.5} * wk.dot(Q_inv * wk);
        }

        // Measurement cost: sum 0.5 * v_k^T R_inv v_k
        for(int k = 0; k <= Ni; ++k)
        {
            Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + x_offset + k * nx);

            Vector<Scalar, NY> vk = state->z_window[static_cast<std::size_t>(k)] - measurement(xk);
            total += Scalar{0.5} * vk.dot(R_inv * vk);
        }

        // L1 slack penalty for path constraints
        if(has_path_slack)
        {
            for(int k = 0; k <= Ni; ++k)
            {
                Eigen::Map<const Vector<Scalar, NC>> sk(z.data() + path_slack_offset + k * nc);
                total += config.path_penalty.dot(sk);
            }
        }

        return total;
    };

    // Gradient: use finite differences (consistent with NMPC pattern)
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> gradient = [cost](std::span<const Scalar> z, std::span<Scalar> grad) { finite_diff_gradient<Scalar>(cost, z, grad); };

    // Constraints
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> constraints = [=](std::span<const Scalar> z, std::span<Scalar> c)
    {
        // Continuity constraints: x_{k+1} - f(x_k, u_k) = 0
        for(int k = 0; k < Ni; ++k)
        {
            Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + x_offset + k * nx);
            Eigen::Map<const Vector<Scalar, NX>> xk1(z.data() + x_offset + (k + 1) * nx);

            Vector<Scalar, NX> x_next = dynamics(xk, state->u_window[static_cast<std::size_t>(k)]);

            for(int i = 0; i < nx; ++i)
            {
                c[static_cast<std::size_t>(continuity_start + k * nx + i)] = xk1[i] - x_next[i];
            }
        }

        // Path constraints: g(x_k) - s_k <= 0 (soft) or g(x_k) <= 0 (hard)
        if constexpr(NC > 0)
        {
            if(config.path_constraint)
            {
                const auto& g = *config.path_constraint;
                for(int k = 0; k <= Ni; ++k)
                {
                    Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + x_offset + k * nx);

                    Vector<Scalar, NC> gk = g(xk);

                    if(has_path_slack)
                    {
                        Eigen::Map<const Vector<Scalar, NC>> sk(z.data() + path_slack_offset + k * nc);
                        gk -= sk;
                    }

                    for(int i = 0; i < nc; ++i)
                    {
                        c[static_cast<std::size_t>(path_con_start + k * nc + i)] = gk[i];
                    }
                }
            }
        }

        // Measurement residual bounds: |z_k,i - h_i(x_k)| <= bound_i
        // Reformulated as:  (z_k,i - h_i(x_k)) - bound_i <= 0
        //               and -(z_k,i - h_i(x_k)) - bound_i <= 0
        if(has_residual)
        {
            const auto& bound = *config.residual_bound;
            for(int k = 0; k <= Ni; ++k)
            {
                Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + x_offset + k * nx);

                Vector<Scalar, NY> residual = state->z_window[static_cast<std::size_t>(k)] - measurement(xk);

                for(int i = 0; i < ny; ++i)
                {
                    std::size_t base = static_cast<std::size_t>(residual_con_start + k * ny * 2 + i * 2);
                    c[base] = residual[i] - bound[i];
                    c[base + 1] = -residual[i] - bound[i];
                }
            }
        }
    };

    // Variable bounds
    Eigen::VectorX<Scalar> x_lower = Eigen::VectorX<Scalar>::Constant(n_vars, -std::numeric_limits<Scalar>::infinity());
    Eigen::VectorX<Scalar> x_upper = Eigen::VectorX<Scalar>::Constant(n_vars, std::numeric_limits<Scalar>::infinity());

    // State box bounds
    if(config.x_min)
    {
        for(int k = 0; k <= Ni; ++k)
        {
            x_lower.segment(x_offset + k * nx, nx) = *config.x_min;
        }
    }
    if(config.x_max)
    {
        for(int k = 0; k <= Ni; ++k)
        {
            x_upper.segment(x_offset + k * nx, nx) = *config.x_max;
        }
    }

    // Slack non-negativity
    if(has_path_slack)
    {
        x_lower.segment(path_slack_offset, n_path_slack).setZero();
    }

    // Constraint bounds
    Eigen::VectorX<Scalar> c_lower = Eigen::VectorX<Scalar>::Zero(n_constraints);
    Eigen::VectorX<Scalar> c_upper = Eigen::VectorX<Scalar>::Zero(n_constraints);

    // Continuity: equality (c_lower = c_upper = 0, already set)

    // Path constraints: c_lower = -inf, c_upper = 0
    for(int i = path_con_start; i < path_con_start + n_path_con; ++i)
    {
        c_lower[i] = -std::numeric_limits<Scalar>::infinity();
    }

    // Residual bound constraints: c_lower = -inf, c_upper = 0
    for(int i = residual_con_start; i < residual_con_start + n_residual_con; ++i)
    {
        c_lower[i] = -std::numeric_limits<Scalar>::infinity();
    }

    return nlp_problem<Scalar>{.n_vars = n_vars,
                               .n_constraints = n_constraints,
                               .cost = std::move(cost),
                               .gradient = std::move(gradient),
                               .constraints = std::move(constraints),
                               .x_lower = std::move(x_lower),
                               .x_upper = std::move(x_upper),
                               .c_lower = std::move(c_lower),
                               .c_upper = std::move(c_upper)};
}

} // namespace detail
} // namespace ctrlpp

#endif
