#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_ARGMIN_BENCH_METRICS_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_ARGMIN_BENCH_METRICS_H

#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/detail/numerical_diff.h"

#include <Eigen/Core>

#include <algorithm>
#include <cstddef>
#include <ostream>
#include <span>
#include <string_view>
#include <vector>

struct quality_metrics
{
    double objective;
    double max_constraint_violation;
    double gradient_norm;
    bool success;
    int iterations;
    double solve_time_ms;
};

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

    // Max constraint violation
    if(problem.n_constraints > 0)
    {
        auto nc = static_cast<std::size_t>(problem.n_constraints);
        std::vector<double> c(nc);
        problem.constraints(
            std::span<const double>{x_buf.data(), n},
            std::span<double>{c.data(), nc});

        double max_viol = 0.0;
        for(std::size_t i = 0; i < nc; ++i)
        {
            double upper_viol = std::max(0.0, c[i] - problem.c_upper[static_cast<Eigen::Index>(i)]);
            double lower_viol = std::max(0.0, problem.c_lower[static_cast<Eigen::Index>(i)] - c[i]);
            max_viol = std::max(max_viol, std::max(upper_viol, lower_viol));
        }
        m.max_constraint_violation = max_viol;
    }
    else
    {
        m.max_constraint_violation = 0.0;
    }

    return m;
}

inline void write_quality_csv_header(std::ostream& os)
{
    os << "system,solver,algorithm,warm_start,nx,horizon,objective,max_violation,gradient_norm,success,iterations,solve_time_ms\n";
}

inline void write_quality_csv_row(std::ostream& os,
                                  std::string_view system,
                                  std::string_view solver,
                                  std::string_view algorithm,
                                  std::string_view warm_start,
                                  int nx,
                                  int horizon,
                                  const quality_metrics& m)
{
    os << system << ',' << solver << ',' << algorithm << ','
       << warm_start << ',' << nx << ',' << horizon << ','
       << m.objective << ',' << m.max_constraint_violation << ','
       << m.gradient_norm << ',' << (m.success ? 1 : 0) << ','
       << m.iterations << ',' << m.solve_time_ms << '\n';
}

#endif
