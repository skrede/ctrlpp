#ifndef HPP_GUARD_BENCHMARKS_PROFILING_NMPC_PENDULUM_PERF_COMMON_H
#define HPP_GUARD_BENCHMARKS_PROFILING_NMPC_PENDULUM_PERF_COMMON_H

#include "bench_construct.h"

#include "nmpc/pendulum.h"

#include "ctrlpp/nmpc.h"

#include <Eigen/Dense>

#include <span>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <iostream>
#include <charconv>
#include <stdexcept>
#include <string_view>

namespace ctrlpp::argmin_perf
{

struct run_options
{
    int             steps   = 500;
    int             horizon = 20;
    Eigen::Vector2d x0{0.6, 0.0};
    std::string     trace_path;
};

struct trace_row
{
    int    step;
    double theta;
    double omega;
    double u;
    int    iterations;
    double solve_time;
    int    status;
    double stage_cost;
};

inline auto parse_int(std::string_view s) -> int
{
    int v = 0;
    auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
    if(ec != std::errc{})
        throw std::invalid_argument(std::string{"not an integer: "} + std::string{s});
    return v;
}

inline auto parse_double(std::string_view s) -> double
{
    return std::stod(std::string{s});
}

inline auto parse_options(int argc, char** argv) -> run_options
{
    run_options opts;
    for(int i = 1; i < argc; ++i)
    {
        std::string_view arg{argv[i]};
        auto take = [&]() -> std::string_view
        {
            if(i + 1 >= argc)
                throw std::invalid_argument(std::string{"missing value for "} + std::string{arg});
            return std::string_view{argv[++i]};
        };

        if(arg == "--steps")        opts.steps      = parse_int(take());
        else if(arg == "--horizon") opts.horizon    = parse_int(take());
        else if(arg == "--theta0")  opts.x0(0)      = parse_double(take());
        else if(arg == "--omega0")  opts.x0(1)      = parse_double(take());
        else if(arg == "--trace")   opts.trace_path = std::string{take()};
        else
            throw std::invalid_argument(std::string{"unknown flag: "} + std::string{arg});
    }
    return opts;
}

template <typename Solver>
auto run_pendulum_closed_loop(const run_options&      opts,
                              std::vector<trace_row>* trace,
                              std::string_view        label)
    -> std::pair<Eigen::Vector2d, double>
{
    auto config = ctrlpp::bench::problems::nmpc::make_pendulum_config(opts.horizon);

    using dynamics_fn = Eigen::Vector2d (*)(const Eigen::Vector2d&,
                                            const Eigen::Matrix<double, 1, 1>&);
    constexpr dynamics_fn pendulum_fn = ctrlpp::bench::problems::nmpc::pendulum_2;
    auto controller = ctrlpp::bench::built_or_exit(
        ctrlpp::nmpc_dynamic<double, 2, 1, Solver, dynamics_fn>::create(pendulum_fn, config), "controller");

    Eigen::Vector2d x          = opts.x0;
    double          total_cost = 0.0;
    int             solved_k   = 0;

    if(trace)
        trace->reserve(static_cast<std::size_t>(opts.steps));

    std::cout << "PERF_REGION_START " << label << " steps=" << opts.steps
              << " horizon=" << opts.horizon << std::endl;

    for(int k = 0; k < opts.steps; ++k)
    {
        auto u = controller.solve(x);
        if(!u.has_value())
            break;

        const double stage_cost = x.squaredNorm() + u->input.squaredNorm() * 0.1;
        total_cost += stage_cost;

        if(trace)
        {
            const auto diag = controller.diagnostics();
            trace->push_back(trace_row{
                .step       = k,
                .theta      = x(0),
                .omega      = x(1),
                .u          = u->input(0),
                .iterations = diag.iterations,
                .solve_time = diag.solve_time,
                .status     = static_cast<int>(diag.status),
                .stage_cost = stage_cost,
            });
        }

        x = pendulum_fn(x, u->input);
        ++solved_k;
    }

    std::cout << "PERF_REGION_END " << label << " solved=" << solved_k
              << " total_cost=" << total_cost << std::endl;

    return {x, total_cost};
}

inline void dump_trace(const std::string&             path,
                       const std::vector<trace_row>&  trace)
{
    std::ofstream f{path};
    f << "step,theta,omega,u,iterations,solve_time_s,status,stage_cost" << std::endl;
    for(const auto& r : trace)
    {
        f << r.step       << ','
          << r.theta      << ','
          << r.omega      << ','
          << r.u          << ','
          << r.iterations << ','
          << r.solve_time << ','
          << r.status     << ','
          << r.stage_cost << std::endl;
    }
}

}

#endif
