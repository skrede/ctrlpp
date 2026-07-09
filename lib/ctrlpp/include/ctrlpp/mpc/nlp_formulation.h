#ifndef HPP_GUARD_CTRLPP_MPC_NLP_FORMULATION_H
#define HPP_GUARD_CTRLPP_MPC_NLP_FORMULATION_H

/// @brief NLP formulation for nonlinear MPC with multiple shooting.
///
/// @cite rawlings2017 -- Rawlings et al., "Model Predictive Control: Theory, Computation, and Design", 2017, Ch. 8 (NMPC, multiple shooting)
/// @cite diehl2002 -- Diehl, Bock, Schloder et al., "Real-Time Optimization and Nonlinear Model Predictive Control of Processes Governed by DAEs", J. Process Control 12(4), 2002 (real-time iteration scheme)

#include "ctrlpp/types.h"

#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/nmpc_config.h"

#include "ctrlpp/detail/numerical_diff.h"

#include <Eigen/Dense>

#include <span>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>
#include <cstddef>
#include <algorithm>
#include <functional>

namespace ctrlpp
{

template <typename Scalar, std::size_t NX, std::size_t NU>
struct nmpc_formulation_state
{
    Vector<Scalar, NX> x0{Vector<Scalar, NX>::Zero()};
    std::vector<Vector<Scalar, NX>> x_ref{};
    Vector<Scalar, NU> u_prev{Vector<Scalar, NU>::Zero()};
};

namespace detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC = 0, std::size_t NTC = 0, typename Dynamics>
auto build_nmpc_problem(const Dynamics& dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config, std::shared_ptr<nmpc_formulation_state<Scalar, NX, NU>> state) -> nlp_problem<Scalar>
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
    auto cost_fn = [=](std::span<const Scalar> z) -> Scalar
    {
        Scalar total{0};

        for(int k = 0; k < N; ++k)
        {
            Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + x_offset + k * nx);
            Eigen::Map<const Vector<Scalar, NU>> uk(z.data() + u_offset + k * nu);

            if(config.stage_cost)
            {
                total += (*config.stage_cost)(xk, uk);
            }
            else
            {
                const auto& x_ref_k = (static_cast<std::size_t>(k) < state->x_ref.size()) ? state->x_ref[static_cast<std::size_t>(k)] : Vector<Scalar, NX>::Zero();
                Vector<Scalar, NX> dx = xk - x_ref_k;
                total += Scalar{0.5} * dx.dot(config.Q * dx) + Scalar{0.5} * uk.dot(config.R * uk);
            }
        }

        // Terminal cost
        Eigen::Map<const Vector<Scalar, NX>> xN(z.data() + x_offset + N * nx);
        if(config.terminal_cost)
        {
            total += (*config.terminal_cost)(xN);
        }
        else
        {
            const auto Qf = config.Qf.value_or(config.Q);
            const auto& x_ref_N = (static_cast<std::size_t>(N) < state->x_ref.size()) ? state->x_ref[static_cast<std::size_t>(N)] : Vector<Scalar, NX>::Zero();
            Vector<Scalar, NX> dx = xN - x_ref_N;
            total += Scalar{0.5} * dx.dot(Qf * dx);
        }

        // L1 penalty on slack variables
        if(has_path_slack)
        {
            for(int k = 0; k < N; ++k)
            {
                Eigen::Map<const Vector<Scalar, NC>> sk(z.data() + path_slack_offset + k * nc);
                total += config.path_penalty.dot(sk);
            }
        }
        if(has_term_slack)
        {
            Eigen::Map<const Vector<Scalar, NTC>> st(z.data() + term_slack_offset);
            total += config.terminal_penalty.dot(st);
        }

        return total;
    };

    std::function<Scalar(std::span<const Scalar>)> cost = cost_fn;

    // Gradient callback.
    //
    // For the default quadratic tracking objective the gradient is analytic:
    // grad = H z + q, where H is the block-diagonal cost Hessian (per-stage Q,
    // terminal Qf, per-input R; the slack block is zero) and q the linear
    // tracking term (-Q x_ref, -Qf x_ref_N). The L1 slack penalty contributes a
    // constant gradient (path_penalty on each path-slack block, terminal_penalty
    // on the terminal-slack block). Evaluating the blocks directly is exactly the
    // H z + q form without materializing H, and allocates nothing per evaluation.
    //
    // When a custom stage_cost or terminal_cost override is supplied the closed
    // form no longer applies, so the finite-difference gradient is retained for
    // that path. Its perturbation scratch is owned by the callback and sized once.
    const bool cost_is_quadratic = !config.stage_cost && !config.terminal_cost;

    std::function<void(std::span<const Scalar>, std::span<Scalar>)> gradient;
    if(cost_is_quadratic)
    {
        gradient = [=](std::span<const Scalar> z, std::span<Scalar> grad)
        {
            const auto Qf = config.Qf.value_or(config.Q);

            for(int k = 0; k < N; ++k)
            {
                Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + x_offset + k * nx);
                Eigen::Map<const Vector<Scalar, NU>> uk(z.data() + u_offset + k * nu);

                const auto& x_ref_k = (static_cast<std::size_t>(k) < state->x_ref.size()) ? state->x_ref[static_cast<std::size_t>(k)] : Vector<Scalar, NX>::Zero();

                Eigen::Map<Vector<Scalar, NX>> gx(grad.data() + x_offset + k * nx);
                gx = config.Q * (xk - x_ref_k);

                Eigen::Map<Vector<Scalar, NU>> gu(grad.data() + u_offset + k * nu);
                gu = config.R * uk;
            }

            // Terminal state block: Qf (x_N - x_ref_N).
            Eigen::Map<const Vector<Scalar, NX>> xN(z.data() + x_offset + N * nx);
            const auto& x_ref_N = (static_cast<std::size_t>(N) < state->x_ref.size()) ? state->x_ref[static_cast<std::size_t>(N)] : Vector<Scalar, NX>::Zero();
            Eigen::Map<Vector<Scalar, NX>> gxN(grad.data() + x_offset + N * nx);
            gxN = Qf * (xN - x_ref_N);

            // L1 slack penalty gradient (constant in the slack variables).
            if(has_path_slack)
            {
                for(int k = 0; k < N; ++k)
                {
                    Eigen::Map<Vector<Scalar, NC>> gs(grad.data() + path_slack_offset + k * nc);
                    gs = config.path_penalty;
                }
            }
            if(has_term_slack)
            {
                Eigen::Map<Vector<Scalar, NTC>> gst(grad.data() + term_slack_offset);
                gst = config.terminal_penalty;
            }
        };
    }
    else
    {
        gradient = [cost, scratch = std::vector<Scalar>(static_cast<std::size_t>(n_vars))](std::span<const Scalar> z, std::span<Scalar> grad) mutable
        { finite_diff_gradient<Scalar>(cost, z, grad, std::span<Scalar>{scratch.data(), scratch.size()}); };
    }

    // Constraint callback
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> constraints = [=](std::span<const Scalar> z, std::span<Scalar> c)
    {
        // Initial state constraint: z[0..nx] - x0 = 0
        for(int i = 0; i < nx; ++i)
        {
            c[static_cast<std::size_t>(eq_start + i)] = z[static_cast<std::size_t>(x_offset + i)] - state->x0[i];
        }

        // Continuity constraints: z[(k+1)*nx..] - f(xk, uk) = 0
        for(int k = 0; k < N; ++k)
        {
            Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + x_offset + k * nx);
            Eigen::Map<const Vector<Scalar, NU>> uk(z.data() + u_offset + k * nu);

            Vector<Scalar, NX> x_next = dynamics(xk, uk);

            for(int i = 0; i < nx; ++i)
            {
                c[static_cast<std::size_t>(eq_start + (k + 1) * nx + i)] = z[static_cast<std::size_t>(x_offset + (k + 1) * nx + i)] - x_next[i];
            }
        }

        // Rate constraints: one-sided formulation
        if(config.du_max)
        {
            const auto& du_max = *config.du_max;

            for(int k = 0; k < N; ++k)
            {
                Eigen::Map<const Vector<Scalar, NU>> uk(z.data() + u_offset + k * nu);

                Vector<Scalar, NU> uk_prev;
                if(k == 0)
                {
                    uk_prev = state->u_prev;
                }
                else
                {
                    uk_prev = Eigen::Map<const Vector<Scalar, NU>>(z.data() + u_offset + (k - 1) * nu);
                }

                for(int j = 0; j < nu; ++j)
                {
                    Scalar du = uk[j] - uk_prev[j];
                    // (uk - uk_prev) - du_max <= 0
                    c[static_cast<std::size_t>(rate_start + k * nu * 2 + j * 2)] = du - du_max[j];
                    // -(uk - uk_prev) - du_max <= 0
                    c[static_cast<std::size_t>(rate_start + k * nu * 2 + j * 2 + 1)] = -du - du_max[j];
                }
            }
        }

        // Path constraints: g(x_k, u_k) - s_k <= 0 (soft) or g(x_k, u_k) <= 0 (hard)
        if(config.path_constraint)
        {
            const auto& g = *config.path_constraint;
            for(int k = 0; k < N; ++k)
            {
                Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + x_offset + k * nx);
                Eigen::Map<const Vector<Scalar, NU>> uk(z.data() + u_offset + k * nu);

                Vector<Scalar, NC> gk = g(xk, uk);

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

        // Terminal constraints: h(x_N) - s_N <= 0 (soft) or h(x_N) <= 0 (hard)
        if(config.terminal_constraint)
        {
            const auto& h = *config.terminal_constraint;
            Eigen::Map<const Vector<Scalar, NX>> xN(z.data() + x_offset + N * nx);

            Vector<Scalar, NTC> hN = h(xN);

            if(has_term_slack)
            {
                Eigen::Map<const Vector<Scalar, NTC>> st(z.data() + term_slack_offset);
                hN -= st;
            }

            for(int i = 0; i < ntc; ++i)
            {
                c[static_cast<std::size_t>(term_con_start + i)] = hN[i];
            }
        }
    };

    // Variable bounds
    Eigen::VectorX<Scalar> x_lower = Eigen::VectorX<Scalar>::Constant(n_vars, -std::numeric_limits<Scalar>::infinity());
    Eigen::VectorX<Scalar> x_upper = Eigen::VectorX<Scalar>::Constant(n_vars, std::numeric_limits<Scalar>::infinity());

    // State bounds (x0 is NOT pinned by variable bounds; equality constraint handles it).
    // The loop starts at k = 1: x0 is pinned solely by the initial-state equality
    // (z[x_offset+i] - x0[i] = 0 above), so also bounding the x0 block here would be
    // redundant and can render the NLP infeasible when x0 sits on or over a bound.
    // Bounds therefore apply to stages 1..N only.
    for(int k = 1; k <= N; ++k)
    {
        if(config.x_min)
        {
            x_lower.segment(x_offset + k * nx, nx) = *config.x_min;
        }
        if(config.x_max)
        {
            x_upper.segment(x_offset + k * nx, nx) = *config.x_max;
        }
    }

    // Input bounds
    for(int k = 0; k < N; ++k)
    {
        if(config.u_min)
        {
            x_lower.segment(u_offset + k * nu, nu) = *config.u_min;
        }
        if(config.u_max)
        {
            x_upper.segment(u_offset + k * nu, nu) = *config.u_max;
        }
    }

    // Slack variable bounds: s >= 0 (lower = 0, upper = +inf, already set)
    if(has_path_slack)
    {
        x_lower.segment(path_slack_offset, n_path_slack).setZero();
    }
    if(has_term_slack)
    {
        x_lower.segment(term_slack_offset, n_term_slack).setZero();
    }

    // Constraint bounds
    Eigen::VectorX<Scalar> c_lower = Eigen::VectorX<Scalar>::Zero(n_constraints);
    Eigen::VectorX<Scalar> c_upper = Eigen::VectorX<Scalar>::Zero(n_constraints);

    // Equality constraints: c_lower = c_upper = 0 (already set)

    // Inequality constraints (rate): c_lower = -inf, c_upper = 0
    for(int i = rate_start; i < rate_start + n_rate; ++i)
    {
        c_lower[i] = -std::numeric_limits<Scalar>::infinity();
    }

    // Path constraints: c_lower = -inf, c_upper = 0
    for(int i = path_con_start; i < path_con_start + n_path_con; ++i)
    {
        c_lower[i] = -std::numeric_limits<Scalar>::infinity();
    }

    // Terminal constraints: c_lower = -inf, c_upper = 0
    for(int i = term_con_start; i < term_con_start + n_term_con; ++i)
    {
        c_lower[i] = -std::numeric_limits<Scalar>::infinity();
    }

    // Constraint Jacobian callback.
    //
    // Scope (bounded by black-box dynamics): the structurally-analytic entries
    // are written exactly, while entries that depend on the caller-supplied
    // dynamics or constraint callables are finite-differenced. The dynamics_model
    // concept exposes only f(x,u) -> x_next with no Jacobian, and the same holds
    // for the path/terminal constraint functors. Concretely:
    //   * analytic (exact): the initial-state identity rows, the continuity I
    //     block on x_{k+1}, the rate +/-1 blocks, and the slack -1 columns;
    //   * finite-difference: the continuity dynamics blocks df/dx, df/du and the
    //     path/terminal constraint gradients dg/dx, dg/du, dh/dx.
    // The full Jacobian is finite-differenced first (central difference, the same
    // eps^(1/3) magnitude-scaled step as detail::finite_diff_*), then the
    // structural entries are overwritten with their exact values so those rows
    // carry no finite-difference noise.
    //
    // Layout: argmin consumes the raw Jacobian as a column-major
    // (n_constraints x n_vars) matrix flattened into the span, so entry (row i,
    // col j) lives at index i + j * n_constraints. The scratch buffers are owned
    // by the callback and sized once, so per-evaluation assembly does not allocate.
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> constraint_jacobian =
        [=, constraints_fn = constraints,
         z_scratch = std::vector<Scalar>(static_cast<std::size_t>(n_vars)),
         c_plus = std::vector<Scalar>(static_cast<std::size_t>(n_constraints)),
         c_minus = std::vector<Scalar>(static_cast<std::size_t>(n_constraints))](std::span<const Scalar> z, std::span<Scalar> jac) mutable
    {
        const int m = n_constraints;
        const int n = n_vars;
        const auto at = [m](std::span<Scalar> J, int i, int j) -> Scalar&
        { return J[static_cast<std::size_t>(i) + static_cast<std::size_t>(j) * static_cast<std::size_t>(m)]; };

        // Finite-difference the full Jacobian (central difference).
        const auto step_scale = std::cbrt(std::numeric_limits<Scalar>::epsilon());
        std::copy(z.begin(), z.end(), z_scratch.begin());

        for(int j = 0; j < n; ++j)
        {
            const auto jz = static_cast<std::size_t>(j);
            const Scalar h_raw = step_scale * std::max(Scalar{1}, std::abs(z[jz]));
            const Scalar temp = z[jz] + h_raw;
            const Scalar h = temp - z[jz];
            const Scalar orig = z_scratch[jz];

            z_scratch[jz] = orig + h;
            constraints_fn(std::span<const Scalar>{z_scratch.data(), z_scratch.size()}, std::span<Scalar>{c_plus.data(), c_plus.size()});

            z_scratch[jz] = orig - h;
            constraints_fn(std::span<const Scalar>{z_scratch.data(), z_scratch.size()}, std::span<Scalar>{c_minus.data(), c_minus.size()});

            for(int i = 0; i < m; ++i)
            {
                at(jac, i, j) = (c_plus[static_cast<std::size_t>(i)] - c_minus[static_cast<std::size_t>(i)]) / (Scalar{2} * h);
            }

            z_scratch[jz] = orig;
        }

        // Overwrite the structurally-analytic entries with their exact values.

        // Initial-state equality: c[i] = z[x_offset+i] - x0[i]  =>  row i is e_i.
        for(int i = 0; i < nx; ++i)
        {
            for(int j = 0; j < n; ++j)
            {
                at(jac, eq_start + i, j) = Scalar{0};
            }
            at(jac, eq_start + i, x_offset + i) = Scalar{1};
        }

        // Continuity equality: c[(k+1)*nx+i] = z[x_offset+(k+1)*nx+i] - f(xk,uk)[i].
        // The identity block on x_{k+1} is exact; df/dx, df/du stay finite-differenced.
        for(int k = 0; k < N; ++k)
        {
            for(int i = 0; i < nx; ++i)
            {
                at(jac, eq_start + (k + 1) * nx + i, x_offset + (k + 1) * nx + i) = Scalar{1};
            }
        }

        // Rate constraints are fully structural (linear in the inputs):
        //   row (k, j, +):  du - du_max   =>  d/duk[j] = +1,  d/du_{k-1}[j] = -1
        //   row (k, j, -): -du - du_max   =>  d/duk[j] = -1,  d/du_{k-1}[j] = +1
        // At k = 0 the previous input is the constant u_prev (no derivative).
        if(config.du_max)
        {
            for(int k = 0; k < N; ++k)
            {
                for(int j = 0; j < nu; ++j)
                {
                    const int row_p = rate_start + k * nu * 2 + j * 2;
                    const int row_m = row_p + 1;

                    for(int col = 0; col < n; ++col)
                    {
                        at(jac, row_p, col) = Scalar{0};
                        at(jac, row_m, col) = Scalar{0};
                    }

                    at(jac, row_p, u_offset + k * nu + j) = Scalar{1};
                    at(jac, row_m, u_offset + k * nu + j) = Scalar{-1};

                    if(k > 0)
                    {
                        at(jac, row_p, u_offset + (k - 1) * nu + j) = Scalar{-1};
                        at(jac, row_m, u_offset + (k - 1) * nu + j) = Scalar{1};
                    }
                }
            }
        }

        // Path constraint slack columns: c = g(xk,uk) - s_path_k  =>  d/ds_path = -1.
        // dg/dx, dg/du stay finite-differenced.
        if(config.path_constraint && has_path_slack)
        {
            for(int k = 0; k < N; ++k)
            {
                for(int i = 0; i < nc; ++i)
                {
                    at(jac, path_con_start + k * nc + i, path_slack_offset + k * nc + i) = Scalar{-1};
                }
            }
        }

        // Terminal constraint slack columns: c = h(xN) - s_term  =>  d/ds_term = -1.
        // dh/dx stays finite-differenced.
        if(config.terminal_constraint && has_term_slack)
        {
            for(int i = 0; i < ntc; ++i)
            {
                at(jac, term_con_start + i, term_slack_offset + i) = Scalar{-1};
            }
        }
    };

    return nlp_problem<Scalar>{.n_vars = n_vars,
                               .n_constraints = n_constraints,
                               .cost = std::move(cost),
                               .gradient = std::move(gradient),
                               .constraints = std::move(constraints),
                               .constraint_jacobian = std::move(constraint_jacobian),
                               .x_lower = std::move(x_lower),
                               .x_upper = std::move(x_upper),
                               .c_lower = std::move(c_lower),
                               .c_upper = std::move(c_upper)};
}

}
}

#endif
