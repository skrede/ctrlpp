#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_ARGMIN_BENCH_METRICS_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_ARGMIN_BENCH_METRICS_H

#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/argmin_policies.h"
#include "ctrlpp/detail/numerical_diff.h"

#include <Eigen/Core>

#include <span>
#include <vector>
#include <ostream>
#include <cstddef>
#include <algorithm>
#include <string_view>
#include <limits>

struct quality_metrics
{
    double objective;
    double max_constraint_violation;
    double gradient_norm;
    bool success;
    int iterations;
    double solve_time_ms;
};

/// How far the decision vector is from satisfying the constraints the problem
/// poses on it, re-evaluated from the problem rather than taken from whatever
/// the solver reported about its own answer.
inline auto max_constraint_violation(const ctrlpp::nlp_problem<double>& problem,
                                     std::span<const double> x) -> double
{
    if(problem.n_constraints <= 0)
        return 0.0;

    auto nc = static_cast<std::size_t>(problem.n_constraints);
    std::vector<double> c(nc);
    problem.constraints(x, std::span<double>{c.data(), nc});

    double worst = 0.0;
    for(std::size_t i = 0; i < nc; ++i)
    {
        double above = c[i] - problem.c_upper[static_cast<Eigen::Index>(i)];
        double below = problem.c_lower[static_cast<Eigen::Index>(i)] - c[i];
        worst = std::max({worst, above, below});
    }
    return worst;
}

inline auto compute_quality_metrics(const ctrlpp::nlp_problem<double>& problem,
                                    const ctrlpp::nlp_result<double>& result) -> quality_metrics
{
    quality_metrics m{};
    m.objective = result.objective;
    m.success = (result.status == ctrlpp::solve_status::optimal);
    m.iterations = result.iterations;
    m.solve_time_ms = result.solve_time * 1000.0;

    // Gradient norm via finite differences
    auto n = static_cast<std::size_t>(problem.n_vars);
    std::vector<double> grad(n);
    std::vector<double> x_buf(result.x.data(), result.x.data() + result.x.size());
    ctrlpp::detail::finite_diff_gradient<double>(
        problem.cost,
        std::span<const double>{x_buf.data(), n},
        std::span<double>{grad.data(), n});
    m.gradient_norm = Eigen::Map<Eigen::VectorXd>(grad.data(), static_cast<Eigen::Index>(n)).norm();

    m.max_constraint_violation = max_constraint_violation(problem, std::span<const double>{x_buf.data(), n});

    return m;
}

/// The own criterion of a predictive-control arm: the trajectory it returned is
/// claimed to satisfy the shooting continuity of the plant, and this is how far
/// from doing so it is. Solvers exiting on a budget leave last_solution() empty,
/// which the same guard compute_gradient_norm carries reports rather than reads
/// past.
template <typename NmpcType>
auto compute_constraint_violation(const NmpcType& controller) -> double
{
    const auto& problem = controller.problem();
    const auto& z = controller.last_solution();
    auto n = static_cast<std::size_t>(problem.n_vars);
    if(static_cast<std::size_t>(z.size()) != n)
        return std::numeric_limits<double>::quiet_NaN();

    return max_constraint_violation(problem, std::span<const double>{z.data(), n});
}

template <typename Scalar, std::size_t NX, std::size_t NU, typename NmpcType>
auto compute_gradient_norm(const NmpcType& controller) -> double
{
    const auto& prob = controller.problem();
    const auto& z = controller.last_solution();
    auto n = static_cast<std::size_t>(prob.n_vars);

    // ctrlpp::nmpc::solve_impl writes m_last_solution only on solve_status::optimal
    // / solved_inaccurate. Solvers exiting time_limit / max_iterations / etc. leave
    // last_solution() at its default-empty state, so a span sized to prob.n_vars
    // would overrun the buffer in finite_diff_gradient. Sentinel out instead.
    if(static_cast<std::size_t>(z.size()) != n)
        return std::numeric_limits<double>::quiet_NaN();

    std::vector<double> grad(n);
    std::vector<double> x_buf(z.data(), z.data() + z.size());
    ctrlpp::detail::finite_diff_gradient<double>(
        prob.cost,
        std::span<const double>{x_buf.data(), n},
        std::span<double>{grad.data(), n});
    return Eigen::Map<Eigen::VectorXd>(grad.data(), static_cast<Eigen::Index>(n)).norm();
}

/// Names the mode for the warm_start column below. A quality record that spells
/// the column itself can state a mode the run did not use, so every writer
/// derives it from the settings the arms were built from.
inline auto warm_start_label(ctrlpp::warm_start_mode mode) -> std::string
{
    switch(mode)
    {
    case ctrlpp::warm_start_mode::cold:        return "cold";
    case ctrlpp::warm_start_mode::primal_only: return "primal_only";
    case ctrlpp::warm_start_mode::curvature:   return "curvature";
    }
    return "unknown";
}

/// The warm_start cell for an arm whose library exposes no warm-start channel.
/// Such an arm runs identically in every pass of a warm-start sweep, so naming
/// the sweep's mode on its row would report a setting it never received; the
/// pass it belongs to is carried by the variant column instead.
constexpr char const* warm_start_absent = "n/a";

/// Whether an arm's family exposes a warm-start channel at all. The family name
/// is what the variant loops already carry, and it is the property that decides
/// the cell: a library with no such channel cannot be told a mode, so reporting
/// one against it would describe the sweep rather than the run.
inline auto has_warm_start_channel(std::string_view family) -> bool
{
    return family == "argmin";
}

inline void write_quality_csv_header(std::ostream& os)
{
    os << "system,solver,algorithm,variant,warm_start,nx,horizon,objective,max_violation,gradient_norm,"
          "success,iterations,solve_time_ms\n";
}

/// The variant cell names the pass this row belongs to; the warm_start cell
/// names what this arm actually did. They are separate because a pass is not a
/// setting: a sweep over warm-start modes, a Jacobian source and a
/// convergence-rate census all discriminate rows without any of them being a
/// warm-start mode, and one column carrying all four cannot be read.
inline void write_quality_csv_row(std::ostream& os,
                                  std::string_view system,
                                  std::string_view solver,
                                  std::string_view algorithm,
                                  std::string_view variant,
                                  std::string_view warm_start,
                                  int nx,
                                  int horizon,
                                  const quality_metrics& m)
{
    os << system << ',' << solver << ',' << algorithm << ','
       << variant << ',' << warm_start << ',' << nx << ',' << horizon << ','
       << m.objective << ',' << m.max_constraint_violation << ','
       << m.gradient_norm << ',' << (m.success ? 1 : 0) << ','
       << m.iterations << ',' << m.solve_time_ms << '\n';
    os.flush();
}

#endif
