#include "bench_metrics.h"

#include "ctrlpp/nmpc.h"
#include "ctrlpp/mpc/nlopt_solver.h"
#include "ctrlpp/mpc/argmin_solver.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <string>

namespace
{

// ---------------------------------------------------------------------------
// Budget CSV helpers (local to this benchmark)
// ---------------------------------------------------------------------------

void write_budget_csv_header(std::ostream& os)
{
    os << "system,solver,algorithm,budget,nx,horizon,objective,max_violation,"
          "gradient_norm,success,iterations,solve_time_ms\n";
}

void write_budget_csv_row(std::ostream& os,
                          std::string_view system,
                          std::string_view solver,
                          std::string_view algorithm,
                          std::string_view budget,
                          int nx,
                          int horizon,
                          const quality_metrics& m)
{
    os << system << ',' << solver << ',' << algorithm << ','
       << budget << ',' << nx << ',' << horizon << ','
       << m.objective << ',' << m.max_constraint_violation << ','
       << m.gradient_norm << ',' << (m.success ? 1 : 0) << ','
       << m.iterations << ',' << m.solve_time_ms << '\n';
}

// ---------------------------------------------------------------------------
// Dynamics definitions
// ---------------------------------------------------------------------------

constexpr double di2_dt = 0.1;
auto double_integrator_2 = [](const Eigen::Vector2d& x,
                              const Eigen::Matrix<double, 1, 1>& u) -> Eigen::Vector2d
{
    return Eigen::Vector2d{x(0) + di2_dt * x(1), x(1) + di2_dt * u(0)};
};

auto pendulum_2 = [](const Eigen::Vector2d& x,
                     const Eigen::Matrix<double, 1, 1>& u) -> Eigen::Vector2d
{
    constexpr double dt = 0.05;
    constexpr double g = 9.81;
    constexpr double l = 1.0;
    double theta = x(0);
    double omega = x(1);
    double alpha = -g / l * std::sin(theta) + u(0);
    return Eigen::Vector2d{theta + dt * omega, omega + dt * alpha};
};

constexpr double di4_dt = 0.1;
auto double_integrator_4 = [](const Eigen::Vector4d& x,
                              const Eigen::Vector2d& u) -> Eigen::Vector4d
{
    return Eigen::Vector4d{
        x(0) + di4_dt * x(1),
        x(1) + di4_dt * u(0),
        x(2) + di4_dt * x(3),
        x(3) + di4_dt * u(1)};
};

constexpr double di8_dt = 0.1;
using Vec8 = Eigen::Matrix<double, 8, 1>;
using Vec4 = Eigen::Vector4d;

auto double_integrator_8 = [](const Vec8& x, const Vec4& u) -> Vec8
{
    Vec8 xn;
    xn(0) = x(0) + di8_dt * x(1);
    xn(1) = x(1) + di8_dt * u(0);
    xn(2) = x(2) + di8_dt * x(3);
    xn(3) = x(3) + di8_dt * u(1);
    xn(4) = x(4) + di8_dt * x(5);
    xn(5) = x(5) + di8_dt * u(2);
    xn(6) = x(6) + di8_dt * x(7);
    xn(7) = x(7) + di8_dt * u(3);
    return xn;
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
// Type aliases
// ---------------------------------------------------------------------------

using NloptSolver = ctrlpp::nlopt_solver<double>;
using ArgminSlsqp = ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp>;

// ---------------------------------------------------------------------------
// Step-budget sweep runner
// ---------------------------------------------------------------------------

template <std::size_t NX, std::size_t NU, typename Dynamics>
void run_step_budget(const std::string& system_name,
                     Dynamics dynamics,
                     int horizon,
                     std::ostream& quality_csv)
{
    auto config = make_nmpc_config<NX, NU>(horizon);

    Eigen::Matrix<double, NX, 1> x0 = Eigen::Matrix<double, NX, 1>::Zero();
    x0(0) = 1.0;

    constexpr std::array budgets = {1, 2, 5, 10, 20, 50, 100};

    // Budget-limited solves via argmin_solver::step()
    for(int budget : budgets)
    {
        // Fresh controller to get a clean nlp_problem, then fresh solver
        ctrlpp::nmpc<double, NX, NU, ArgminSlsqp, Dynamics> controller{dynamics, config};
        const auto& problem = controller.problem();

        ArgminSlsqp solver{};
        solver.setup(problem);

        int n_vars = problem.n_vars;
        ctrlpp::nlp_update<double> update{};
        update.x0 = Eigen::VectorXd::Zero(n_vars);

        auto result = solver.step(update, budget);

        auto qm = compute_quality_metrics(problem, result);

        write_budget_csv_row(quality_csv, system_name, "argmin", "slsqp",
                             std::to_string(budget),
                             static_cast<int>(NX), horizon, qm);
    }

    // Full argmin solve
    {
        ctrlpp::nmpc<double, NX, NU, ArgminSlsqp, Dynamics> controller{dynamics, config};
        const auto& problem = controller.problem();

        ArgminSlsqp solver{};
        solver.setup(problem);

        int n_vars = problem.n_vars;
        ctrlpp::nlp_update<double> update{};
        update.x0 = Eigen::VectorXd::Zero(n_vars);

        auto result = solver.solve(update);

        auto qm = compute_quality_metrics(problem, result);

        write_budget_csv_row(quality_csv, system_name, "argmin", "slsqp",
                             "full",
                             static_cast<int>(NX), horizon, qm);
    }

    // NLopt full-solve baseline
    {
        ctrlpp::nmpc<double, NX, NU, NloptSolver, Dynamics> controller{dynamics, config};
        controller.solve(x0);
        auto diag = controller.diagnostics();
        auto grad = compute_gradient_norm<double, NX, NU>(controller);

        quality_metrics qm{
            .objective = diag.cost,
            .max_constraint_violation = diag.max_constraint_violation,
            .gradient_norm = grad,
            .success = (diag.status == ctrlpp::solve_status::optimal),
            .iterations = diag.iterations,
            .solve_time_ms = diag.solve_time * 1000.0,
        };

        write_budget_csv_row(quality_csv, system_name, "nlopt", "slsqp",
                             "full_nlopt",
                             static_cast<int>(NX), horizon, qm);
    }
}

}

int main()
{
    std::ofstream quality_csv("bench_step_budget_quality.csv");
    write_budget_csv_header(quality_csv);

    // Problem instances
    run_step_budget<2, 1>("double_integrator", double_integrator_2, 10, quality_csv);
    run_step_budget<4, 2>("double_integrator", double_integrator_4, 10, quality_csv);
    run_step_budget<8, 4>("double_integrator", double_integrator_8, 10, quality_csv);
    run_step_budget<2, 1>("pendulum", pendulum_2, 5, quality_csv);
    run_step_budget<4, 2>("double_integrator", double_integrator_4, 20, quality_csv);
    run_step_budget<4, 2>("double_integrator", double_integrator_4, 30, quality_csv);
}
