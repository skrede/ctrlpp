#ifndef HPP_GUARD_CTRLPP_MPC_NLP_FORMULATION_H
#define HPP_GUARD_CTRLPP_MPC_NLP_FORMULATION_H

/// @brief NLP formulation for nonlinear MPC with multiple shooting.
///
/// @cite rawlings2017 -- Rawlings et al., "Model Predictive Control: Theory, Computation, and Design", 2017

#include "ctrlpp/detail/numerical_diff.h"
#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/nmpc_config.h"
#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <span>
#include <vector>

namespace ctrlpp {

template<typename Scalar, std::size_t NX, std::size_t NU>
struct nmpc_formulation_state {
    Vector<Scalar, NX> x0{Vector<Scalar, NX>::Zero()};
    std::vector<Vector<Scalar, NX>> x_ref{};
    Vector<Scalar, NU> u_prev{Vector<Scalar, NU>::Zero()};
};

namespace detail {

template<typename Scalar, std::size_t NX, std::size_t NU,
         std::size_t NC = 0, std::size_t NTC = 0, typename Dynamics>
auto build_nmpc_problem(const Dynamics& dynamics,
                        const nmpc_config<Scalar, NX, NU, NC, NTC>& config,
                        std::shared_ptr<nmpc_formulation_state<Scalar, NX, NU>> state)
    -> nlp_problem<Scalar>
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);
    constexpr int nc = static_cast<int>(NC);
    constexpr int ntc = static_cast<int>(NTC);
    const int N = config.horizon;

    // Decision vector layout:
    //   [x_0, ..., x_N, u_0, ..., u_{N-1}, s_path_0, ..., s_path_{N-1}, s_term]
    // where s_path_k is Vector<NC> slack at node k, s_term is Vector<NTC> slack at terminal.
    // Slack variables only present when soft_constraints=true and NC>0 / NTC>0.
    const bool has_path_slack = config.soft_constraints && NC > 0 && config.path_constraint;
    const bool has_term_slack = config.soft_constraints && NTC > 0 && config.terminal_constraint;
    const int n_path_slack = has_path_slack ? N * nc : 0;
    const int n_term_slack = has_term_slack ? ntc : 0;
    const int n_vars = (N + 1) * nx + N * nu + n_path_slack + n_term_slack;

    // Offsets into decision vector
    const int x_offset = 0;
    const int u_offset = (N + 1) * nx;
    const int path_slack_offset = u_offset + N * nu;
    const int term_slack_offset = path_slack_offset + n_path_slack;

    // Constraint count
    //   Equality: (N+1)*NX (initial state + continuity)
    //   Rate inequality: du_max ? N*NU*2 : 0
    //   Path constraint inequality: path_constraint ? N*NC : 0
    //   Terminal constraint inequality: terminal_constraint ? NTC : 0
    const int n_eq = (N + 1) * nx;
    const int n_rate = config.du_max ? N * nu * 2 : 0;
    const int n_path_con = config.path_constraint ? N * nc : 0;
    const int n_term_con = config.terminal_constraint ? ntc : 0;
    const int n_constraints = n_eq + n_rate + n_path_con + n_term_con;

    // Constraint offset map:
    //   [0, n_eq): equality constraints
    //   [n_eq, n_eq + n_rate): rate inequality
    //   [n_eq + n_rate, n_eq + n_rate + n_path_con): path constraints
    //   [n_eq + n_rate + n_path_con, ...): terminal constraints
    const int eq_start = 0;
    const int rate_start = n_eq;
    const int path_con_start = rate_start + n_rate;
    const int term_con_start = path_con_start + n_path_con;

    // Cost callback
    auto cost_fn = [=](std::span<const Scalar> z) -> Scalar {
        Scalar total{0};

        for (int k = 0; k < N; ++k) {
            Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + x_offset + k * nx);
            Eigen::Map<const Vector<Scalar, NU>> uk(z.data() + u_offset + k * nu);

            if (config.stage_cost) {
                total += (*config.stage_cost)(xk, uk);
            } else {
                const auto& x_ref_k = (static_cast<std::size_t>(k) < state->x_ref.size())
                    ? state->x_ref[static_cast<std::size_t>(k)]
                    : Vector<Scalar, NX>::Zero();
                Vector<Scalar, NX> dx = xk - x_ref_k;
                total += Scalar{0.5} * dx.dot(config.Q * dx)
                       + Scalar{0.5} * uk.dot(config.R * uk);
            }
        }

        // Terminal cost
        Eigen::Map<const Vector<Scalar, NX>> xN(z.data() + x_offset + N * nx);
        if (config.terminal_cost) {
            total += (*config.terminal_cost)(xN);
        } else {
            const auto Qf = config.Qf.value_or(config.Q);
            const auto& x_ref_N = (static_cast<std::size_t>(N) < state->x_ref.size())
                ? state->x_ref[static_cast<std::size_t>(N)]
                : Vector<Scalar, NX>::Zero();
            Vector<Scalar, NX> dx = xN - x_ref_N;
            total += Scalar{0.5} * dx.dot(Qf * dx);
        }

        // L1 penalty on slack variables
        if (has_path_slack) {
            for (int k = 0; k < N; ++k) {
                Eigen::Map<const Vector<Scalar, NC>> sk(
                    z.data() + path_slack_offset + k * nc);
                total += config.path_penalty.dot(sk);
            }
        }
        if (has_term_slack) {
            Eigen::Map<const Vector<Scalar, NTC>> st(z.data() + term_slack_offset);
            total += config.terminal_penalty.dot(st);
        }

        return total;
    };

    std::function<Scalar(std::span<const Scalar>)> cost = cost_fn;

    // Gradient callback via finite differences
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> gradient =
        [cost](std::span<const Scalar> z, std::span<Scalar> grad) {
            finite_diff_gradient<Scalar>(cost, z, grad);
        };

    // Constraint callback
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> constraints =
        [=](std::span<const Scalar> z, std::span<Scalar> c) {
            // Initial state constraint: z[0..nx] - x0 = 0
            for (int i = 0; i < nx; ++i) {
                c[static_cast<std::size_t>(eq_start + i)] =
                    z[static_cast<std::size_t>(x_offset + i)] - state->x0[i];
            }

            // Continuity constraints: z[(k+1)*nx..] - f(xk, uk) = 0
            for (int k = 0; k < N; ++k) {
                Eigen::Map<const Vector<Scalar, NX>> xk(
                    z.data() + x_offset + k * nx);
                Eigen::Map<const Vector<Scalar, NU>> uk(
                    z.data() + u_offset + k * nu);

                Vector<Scalar, NX> x_next = dynamics(xk, uk);

                for (int i = 0; i < nx; ++i) {
                    c[static_cast<std::size_t>(eq_start + (k + 1) * nx + i)] =
                        z[static_cast<std::size_t>(x_offset + (k + 1) * nx + i)] - x_next[i];
                }
            }

            // Rate constraints: one-sided formulation
            if (config.du_max) {
                const auto& du_max = *config.du_max;

                for (int k = 0; k < N; ++k) {
                    Eigen::Map<const Vector<Scalar, NU>> uk(
                        z.data() + u_offset + k * nu);

                    Vector<Scalar, NU> uk_prev;
                    if (k == 0) {
                        uk_prev = state->u_prev;
                    } else {
                        uk_prev = Eigen::Map<const Vector<Scalar, NU>>(
                            z.data() + u_offset + (k - 1) * nu);
                    }

                    for (int j = 0; j < nu; ++j) {
                        Scalar du = uk[j] - uk_prev[j];
                        // (uk - uk_prev) - du_max <= 0
                        c[static_cast<std::size_t>(rate_start + k * nu * 2 + j * 2)] =
                            du - du_max[j];
                        // -(uk - uk_prev) - du_max <= 0
                        c[static_cast<std::size_t>(rate_start + k * nu * 2 + j * 2 + 1)] =
                            -du - du_max[j];
                    }
                }
            }

            // Path constraints: g(x_k, u_k) - s_k <= 0 (soft) or g(x_k, u_k) <= 0 (hard)
            if (config.path_constraint) {
                const auto& g = *config.path_constraint;
                for (int k = 0; k < N; ++k) {
                    Eigen::Map<const Vector<Scalar, NX>> xk(
                        z.data() + x_offset + k * nx);
                    Eigen::Map<const Vector<Scalar, NU>> uk(
                        z.data() + u_offset + k * nu);

                    Vector<Scalar, NC> gk = g(xk, uk);

                    if (has_path_slack) {
                        Eigen::Map<const Vector<Scalar, NC>> sk(
                            z.data() + path_slack_offset + k * nc);
                        gk -= sk;
                    }

                    for (int i = 0; i < nc; ++i) {
                        c[static_cast<std::size_t>(path_con_start + k * nc + i)] = gk[i];
                    }
                }
            }

            // Terminal constraints: h(x_N) - s_N <= 0 (soft) or h(x_N) <= 0 (hard)
            if (config.terminal_constraint) {
                const auto& h = *config.terminal_constraint;
                Eigen::Map<const Vector<Scalar, NX>> xN(
                    z.data() + x_offset + N * nx);

                Vector<Scalar, NTC> hN = h(xN);

                if (has_term_slack) {
                    Eigen::Map<const Vector<Scalar, NTC>> st(
                        z.data() + term_slack_offset);
                    hN -= st;
                }

                for (int i = 0; i < ntc; ++i) {
                    c[static_cast<std::size_t>(term_con_start + i)] = hN[i];
                }
            }
        };

    // Variable bounds
    Eigen::VectorX<Scalar> x_lower = Eigen::VectorX<Scalar>::Constant(
        n_vars, -std::numeric_limits<Scalar>::infinity());
    Eigen::VectorX<Scalar> x_upper = Eigen::VectorX<Scalar>::Constant(
        n_vars, std::numeric_limits<Scalar>::infinity());

    // State bounds (x0 is NOT pinned by variable bounds; equality constraint handles it)
    for (int k = 0; k <= N; ++k) {
        if (config.x_min) {
            x_lower.segment(x_offset + k * nx, nx) = *config.x_min;
        }
        if (config.x_max) {
            x_upper.segment(x_offset + k * nx, nx) = *config.x_max;
        }
    }

    // Input bounds
    for (int k = 0; k < N; ++k) {
        if (config.u_min) {
            x_lower.segment(u_offset + k * nu, nu) = *config.u_min;
        }
        if (config.u_max) {
            x_upper.segment(u_offset + k * nu, nu) = *config.u_max;
        }
    }

    // Slack variable bounds: s >= 0 (lower = 0, upper = +inf, already set)
    if (has_path_slack) {
        x_lower.segment(path_slack_offset, n_path_slack).setZero();
    }
    if (has_term_slack) {
        x_lower.segment(term_slack_offset, n_term_slack).setZero();
    }

    // Constraint bounds
    Eigen::VectorX<Scalar> c_lower = Eigen::VectorX<Scalar>::Zero(n_constraints);
    Eigen::VectorX<Scalar> c_upper = Eigen::VectorX<Scalar>::Zero(n_constraints);

    // Equality constraints: c_lower = c_upper = 0 (already set)

    // Inequality constraints (rate): c_lower = -inf, c_upper = 0
    for (int i = rate_start; i < rate_start + n_rate; ++i) {
        c_lower[i] = -std::numeric_limits<Scalar>::infinity();
    }

    // Path constraints: c_lower = -inf, c_upper = 0
    for (int i = path_con_start; i < path_con_start + n_path_con; ++i) {
        c_lower[i] = -std::numeric_limits<Scalar>::infinity();
    }

    // Terminal constraints: c_lower = -inf, c_upper = 0
    for (int i = term_con_start; i < term_con_start + n_term_con; ++i) {
        c_lower[i] = -std::numeric_limits<Scalar>::infinity();
    }

    return nlp_problem<Scalar>{
        .n_vars = n_vars,
        .n_constraints = n_constraints,
        .cost = std::move(cost),
        .gradient = std::move(gradient),
        .constraints = std::move(constraints),
        .x_lower = std::move(x_lower),
        .x_upper = std::move(x_upper),
        .c_lower = std::move(c_lower),
        .c_upper = std::move(c_upper)
    };
}

}
}

#endif
