#include "bench_metrics.h"

#include "ctrlpp/nmpc.h"
#include "ctrlpp/types.h"
#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/nmpc_config.h"
#include "ctrlpp/mpc/argmin_solver.h"
#include "ctrlpp/mpc/nlp_formulation.h"
#include "ctrlpp/detail/numerical_diff.h"

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace
{

constexpr char const* comma_csv_tpl = R"TEMPLATE(
"title","name","unit","batch","elapsed","error%","instructions","branches","branch_misses","total"
{{#result}}"{{title}}","{{name}}","{{unit}}",{{batch}},{{median(elapsed)}},{{medianAbsolutePercentError(elapsed)}},{{median(instructions)}},{{median(branchinstructions)}},{{median(branchmisses)}},{{sumProduct(iterations, elapsed)}}
{{/result}})TEMPLATE";

// ---------------------------------------------------------------------------
// Differentiable dynamics definitions
// ---------------------------------------------------------------------------

struct differentiable_double_integrator_2
{
    static constexpr double dt = 0.1;

    auto operator()(const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>& u) const -> Eigen::Vector2d
    {
        return {x(0) + dt * x(1), x(1) + dt * u(0)};
    }

    auto jacobian_x(const Eigen::Vector2d&, const Eigen::Matrix<double, 1, 1>&) const -> Eigen::Matrix2d
    {
        Eigen::Matrix2d A;
        A << 1.0, dt, 0.0, 1.0;
        return A;
    }

    auto jacobian_u(const Eigen::Vector2d&, const Eigen::Matrix<double, 1, 1>&) const -> Eigen::Matrix<double, 2, 1>
    {
        return Eigen::Matrix<double, 2, 1>{0.0, dt};
    }
};

struct differentiable_pendulum_2
{
    static constexpr double dt = 0.05;
    static constexpr double g = 9.81;
    static constexpr double l = 1.0;

    auto operator()(const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>& u) const -> Eigen::Vector2d
    {
        double theta = x(0);
        double omega = x(1);
        double alpha = -g / l * std::sin(theta) + u(0);
        return {theta + dt * omega, omega + dt * alpha};
    }

    auto jacobian_x(const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>&) const -> Eigen::Matrix2d
    {
        Eigen::Matrix2d A;
        A << 1.0, dt, -g / l * std::cos(x(0)) * dt, 1.0;
        return A;
    }

    auto jacobian_u(const Eigen::Vector2d&, const Eigen::Matrix<double, 1, 1>&) const -> Eigen::Matrix<double, 2, 1>
    {
        return Eigen::Matrix<double, 2, 1>{0.0, dt};
    }
};

// ---------------------------------------------------------------------------
// NMPC config factory
// ---------------------------------------------------------------------------

template <std::size_t NX, std::size_t NU>
auto make_nmpc_config(int horizon) -> ctrlpp::nmpc_config<double, NX, NU>
{
    return {
        .horizon = horizon,
        .Q = Eigen::Matrix<double, NX, NX>::Identity(),
        .R = Eigen::Matrix<double, NU, NU>::Identity() * 0.1,
    };
}

// ---------------------------------------------------------------------------
// Type alias
// ---------------------------------------------------------------------------

using ArgminSlsqp = ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp>;

// ---------------------------------------------------------------------------
// Analytic constraint Jacobian builder
// ---------------------------------------------------------------------------

template <std::size_t NX, std::size_t NU, typename Dynamics>
auto build_analytic_constraint_jacobian(
    const Dynamics& dynamics,
    int horizon,
    std::shared_ptr<ctrlpp::nmpc_formulation_state<double, NX, NU>> state)
    -> std::function<void(std::span<const double>, std::span<double>)>
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);

    return [=](std::span<const double> z, std::span<double> jac_flat)
    {
        const int N = horizon;
        const int n_vars = (N + 1) * nx + N * nu;
        const int n_constraints = (N + 1) * nx;
        const int x_offset = 0;
        const int u_offset = (N + 1) * nx;

        std::fill(jac_flat.begin(), jac_flat.end(), 0.0);

        // Initial state constraint: x_0 - x0_param = 0
        // Jacobian rows 0..NX-1: identity block at x_0 columns
        for(int i = 0; i < nx; ++i)
        {
            jac_flat[static_cast<std::size_t>(i * n_vars + x_offset + i)] = 1.0;
        }

        // Continuity constraints at node k: x_{k+1} - f(x_k, u_k) = 0
        for(int k = 0; k < N; ++k)
        {
            Eigen::Map<const ctrlpp::Vector<double, NX>> xk(z.data() + x_offset + k * nx);
            Eigen::Map<const ctrlpp::Vector<double, NU>> uk(z.data() + u_offset + k * nu);

            auto A = dynamics.jacobian_x(xk, uk);
            auto B = dynamics.jacobian_u(xk, uk);

            int row_base = (k + 1) * nx;

            // d/dx_k: -A_k
            for(int i = 0; i < nx; ++i)
            {
                for(int j = 0; j < nx; ++j)
                {
                    jac_flat[static_cast<std::size_t>((row_base + i) * n_vars + x_offset + k * nx + j)] = -A(i, j);
                }
            }

            // d/dx_{k+1}: I
            for(int i = 0; i < nx; ++i)
            {
                jac_flat[static_cast<std::size_t>((row_base + i) * n_vars + x_offset + (k + 1) * nx + i)] = 1.0;
            }

            // d/du_k: -B_k
            for(int i = 0; i < nx; ++i)
            {
                for(int j = 0; j < nu; ++j)
                {
                    jac_flat[static_cast<std::size_t>((row_base + i) * n_vars + u_offset + k * nu + j)] = -B(i, j);
                }
            }
        }
    };
}

// ---------------------------------------------------------------------------
// Benchmark runner
// ---------------------------------------------------------------------------

template <std::size_t NX, std::size_t NU, typename Dynamics>
void run_jacobian_benchmark(const std::string& system_name,
                            const Dynamics& dynamics,
                            int horizon,
                            ankerl::nanobench::Bench& bench,
                            std::ostream& quality_csv)
{
    auto config = make_nmpc_config<NX, NU>(horizon);
    auto state = std::make_shared<ctrlpp::nmpc_formulation_state<double, NX, NU>>();
    state->x_ref.resize(static_cast<std::size_t>(horizon + 1), ctrlpp::Vector<double, NX>::Zero());

    Eigen::Matrix<double, NX, 1> x0 = Eigen::Matrix<double, NX, 1>::Zero();
    x0(0) = 1.0;
    state->x0 = x0;

    auto title = system_name + " NX=" + std::to_string(NX)
               + " N=" + std::to_string(horizon);

    // Build the NLP problem (FD gradient, no analytic Jacobian yet)
    auto problem_fd = ctrlpp::detail::build_nmpc_problem<double, NX, NU>(dynamics, config, state);

    // Build analytic Jacobian variant
    auto problem_analytic = problem_fd;
    problem_analytic.constraint_jacobian = build_analytic_constraint_jacobian<NX, NU>(dynamics, horizon, state);

    // FD timing
    {
        ctrlpp::argmin_settings<double> cfg{};
        ArgminSlsqp solver{cfg};
        solver.setup(problem_fd);

        bench.warmup(50).minEpochIterations(50).title(title)
            .run("fd_jacobian",
                 [&]
                 {
                     ctrlpp::nlp_update<double> upd{.x0 = Eigen::VectorXd::Zero(problem_fd.n_vars)};
                     auto result = solver.solve(upd);
                     ankerl::nanobench::doNotOptimizeAway(result);
                 });
    }

    // Analytic timing
    {
        ctrlpp::argmin_settings<double> cfg{};
        ArgminSlsqp solver{cfg};
        solver.setup(problem_analytic);

        bench.run("analytic_jacobian",
                  [&]
                  {
                      ctrlpp::nlp_update<double> upd{.x0 = Eigen::VectorXd::Zero(problem_analytic.n_vars)};
                      auto result = solver.solve(upd);
                      ankerl::nanobench::doNotOptimizeAway(result);
                  });
    }

    // Quality: single solve for metrics
    {
        ctrlpp::argmin_settings<double> cfg{};
        ArgminSlsqp solver_fd{cfg};
        solver_fd.setup(problem_fd);
        ctrlpp::nlp_update<double> upd{.x0 = Eigen::VectorXd::Zero(problem_fd.n_vars)};
        auto result_fd = solver_fd.solve(upd);
        auto qm_fd = compute_quality_metrics(problem_fd, result_fd);

        ArgminSlsqp solver_an{cfg};
        solver_an.setup(problem_analytic);
        auto result_an = solver_an.solve(upd);
        auto qm_an = compute_quality_metrics(problem_analytic, result_an);

        write_quality_csv_row(quality_csv, system_name, "argmin", "slsqp", "fd_jacobian",
                              static_cast<int>(NX), horizon, qm_fd);
        write_quality_csv_row(quality_csv, system_name, "argmin", "slsqp", "analytic_jacobian",
                              static_cast<int>(NX), horizon, qm_an);
    }
}

}

int main()
{
    ankerl::nanobench::Bench bench;
    bench.performanceCounters(true).relative(true);

    std::ofstream timing_csv("bench_jacobian_timing.csv");
    std::ofstream quality_csv("bench_jacobian_quality.csv");
    write_quality_csv_header(quality_csv);

    differentiable_double_integrator_2 di2;
    differentiable_pendulum_2 pend2;

    run_jacobian_benchmark<2, 1>("double_integrator", di2, 10, bench, quality_csv);
    run_jacobian_benchmark<2, 1>("double_integrator", di2, 20, bench, quality_csv);
    run_jacobian_benchmark<2, 1>("pendulum", pend2, 10, bench, quality_csv);
    run_jacobian_benchmark<2, 1>("pendulum", pend2, 5, bench, quality_csv);

    bench.render(comma_csv_tpl, timing_csv);
}
