#ifndef HPP_GUARD_CTRLPP_MPC_NLP_FORMULATION_H
#define HPP_GUARD_CTRLPP_MPC_NLP_FORMULATION_H

#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/nmpc_config.h"
#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
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

template<typename Scalar>
void finite_diff_gradient(const std::function<Scalar(std::span<const Scalar>)>& f,
                          std::span<const Scalar> z,
                          std::span<Scalar> grad)
{
    const auto eps = std::sqrt(std::numeric_limits<Scalar>::epsilon());
    const auto n = z.size();

    std::vector<Scalar> z_mut(z.begin(), z.end());

    for (std::size_t j = 0; j < n; ++j) {
        const Scalar h = eps * std::max(Scalar{1}, std::abs(z[j]));
        const Scalar orig = z_mut[j];

        z_mut[j] = orig + h;
        const Scalar f_plus = f(std::span<const Scalar>{z_mut.data(), n});

        z_mut[j] = orig - h;
        const Scalar f_minus = f(std::span<const Scalar>{z_mut.data(), n});

        grad[j] = (f_plus - f_minus) / (Scalar{2} * h);
        z_mut[j] = orig;
    }
}

template<typename Scalar>
void finite_diff_jacobian(const std::function<void(std::span<const Scalar>, std::span<Scalar>)>& c,
                          int n_constraints,
                          std::span<const Scalar> z,
                          std::span<Scalar> jac)
{
    const auto eps = std::sqrt(std::numeric_limits<Scalar>::epsilon());
    const auto n = z.size();
    const auto m = static_cast<std::size_t>(n_constraints);

    std::vector<Scalar> z_mut(z.begin(), z.end());
    std::vector<Scalar> c_plus(m);
    std::vector<Scalar> c_minus(m);

    for (std::size_t j = 0; j < n; ++j) {
        const Scalar h = eps * std::max(Scalar{1}, std::abs(z[j]));
        const Scalar orig = z_mut[j];

        z_mut[j] = orig + h;
        c(std::span<const Scalar>{z_mut.data(), n},
          std::span<Scalar>{c_plus.data(), m});

        z_mut[j] = orig - h;
        c(std::span<const Scalar>{z_mut.data(), n},
          std::span<Scalar>{c_minus.data(), m});

        for (std::size_t i = 0; i < m; ++i) {
            jac[i * n + j] = (c_plus[i] - c_minus[i]) / (Scalar{2} * h);
        }

        z_mut[j] = orig;
    }
}

template<typename Scalar, std::size_t NX, std::size_t NU, typename Dynamics>
auto build_nmpc_problem(const Dynamics& dynamics,
                        const nmpc_config<Scalar, NX, NU>& config,
                        std::shared_ptr<nmpc_formulation_state<Scalar, NX, NU>> state)
    -> nlp_problem<Scalar>
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);
    const int N = config.horizon;
    const int n_vars = (N + 1) * nx + N * nu;

    // Constraint count
    const int n_eq = (N + 1) * nx;
    const int n_ineq = config.du_max ? N * nu * 2 : 0;
    const int n_constraints = n_eq + n_ineq;

    // Cost callback
    auto cost_fn = [dynamics, config, state, N, n_vars](std::span<const Scalar> z) -> Scalar {
        Scalar total{0};

        for (int k = 0; k < N; ++k) {
            Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + k * nx);
            Eigen::Map<const Vector<Scalar, NU>> uk(z.data() + (N + 1) * nx + k * nu);

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
        Eigen::Map<const Vector<Scalar, NX>> xN(z.data() + N * nx);
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
        [dynamics, config, state, N](std::span<const Scalar> z, std::span<Scalar> c) {
            // Initial state constraint: z[0..nx] - x0 = 0
            for (int i = 0; i < nx; ++i) {
                c[static_cast<std::size_t>(i)] = z[static_cast<std::size_t>(i)] - state->x0[i];
            }

            // Continuity constraints: z[(k+1)*nx..] - f(xk, uk) = 0
            for (int k = 0; k < N; ++k) {
                Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + k * nx);
                Eigen::Map<const Vector<Scalar, NU>> uk(z.data() + (N + 1) * nx + k * nu);

                Vector<Scalar, NX> x_next = dynamics(xk, uk);

                for (int i = 0; i < nx; ++i) {
                    c[static_cast<std::size_t>((k + 1) * nx + i)] =
                        z[static_cast<std::size_t>((k + 1) * nx + i)] - x_next[i];
                }
            }

            // Rate constraints: one-sided formulation
            if (config.du_max) {
                const int n_eq = (N + 1) * nx;
                const auto& du_max = *config.du_max;

                for (int k = 0; k < N; ++k) {
                    Eigen::Map<const Vector<Scalar, NU>> uk(
                        z.data() + (N + 1) * nx + k * nu);

                    Vector<Scalar, NU> uk_prev;
                    if (k == 0) {
                        uk_prev = state->u_prev;
                    } else {
                        uk_prev = Eigen::Map<const Vector<Scalar, NU>>(
                            z.data() + (N + 1) * nx + (k - 1) * nu);
                    }

                    for (int j = 0; j < nu; ++j) {
                        Scalar du = uk[j] - uk_prev[j];
                        // (uk - uk_prev) - du_max <= 0
                        c[static_cast<std::size_t>(n_eq + k * nu * 2 + j * 2)] =
                            du - du_max[j];
                        // -(uk - uk_prev) - du_max <= 0
                        c[static_cast<std::size_t>(n_eq + k * nu * 2 + j * 2 + 1)] =
                            -du - du_max[j];
                    }
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
            x_lower.segment(k * nx, nx) = *config.x_min;
        }
        if (config.x_max) {
            x_upper.segment(k * nx, nx) = *config.x_max;
        }
    }

    // Input bounds
    for (int k = 0; k < N; ++k) {
        if (config.u_min) {
            x_lower.segment((N + 1) * nx + k * nu, nu) = *config.u_min;
        }
        if (config.u_max) {
            x_upper.segment((N + 1) * nx + k * nu, nu) = *config.u_max;
        }
    }

    // Constraint bounds
    Eigen::VectorX<Scalar> c_lower = Eigen::VectorX<Scalar>::Zero(n_constraints);
    Eigen::VectorX<Scalar> c_upper = Eigen::VectorX<Scalar>::Zero(n_constraints);

    // Equality constraints: c_lower = c_upper = 0 (already set)

    // Inequality constraints (rate): c_lower = -inf, c_upper = 0
    for (int i = n_eq; i < n_constraints; ++i) {
        c_lower[i] = -std::numeric_limits<Scalar>::infinity();
        // c_upper[i] already 0
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
