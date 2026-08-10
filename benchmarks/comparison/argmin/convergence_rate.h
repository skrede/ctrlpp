#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_ARGMIN_CONVERGENCE_RATE_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_ARGMIN_CONVERGENCE_RATE_H

#include "bench_metrics.h"
#include "bench_construct.h"

#include "ctrlpp/nmpc.h"

#include <Eigen/Core>

#include <string>
#include <vector>
#include <random>
#include <ostream>
#include <cstddef>
#include <utility>

namespace ctrlpp::bench::argmin_arms
{

constexpr int convergence_trials = 100;

/// Every variant is offered the same drawn initial state before the next one is
/// drawn, so the rates below compare the variants rather than their samples.
template <std::size_t NX, std::size_t NU, typename Dynamics, typename Variants>
auto convergence_successes(Dynamics dynamics, const ctrlpp::nmpc_config<double, NX, NU>& config,
                           Variants&& variants) -> std::vector<int>
{
    std::vector<int> successes;
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> draw(-2.0, 2.0);

    for(int trial = 0; trial < convergence_trials; ++trial)
    {
        Eigen::Matrix<double, NX, 1> x0;
        for(std::size_t i = 0; i < NX; ++i)
            x0(static_cast<Eigen::Index>(i)) = draw(rng);
        std::size_t index = 0;
        variants(
            [&](auto solver, char const*, char const*, bool, bool counted)
            {
                const std::size_t slot = index++;
                if(successes.size() <= slot)
                    successes.push_back(0);
                if(!counted)
                    return;
                auto controller = ctrlpp::bench::built_or_exit(
                    ctrlpp::nmpc_dynamic<double, NX, NU, decltype(solver), Dynamics>::create(
                        dynamics, config, std::move(solver)),
                    "arm");
                if(controller.solve(x0).has_value())
                    ++successes[slot];
            });
    }
    return successes;
}

template <std::size_t NX, std::size_t NU, typename Dynamics, typename Variants>
void write_convergence_rates(const std::string& system_name, Dynamics dynamics,
                             const ctrlpp::nmpc_config<double, NX, NU>& config, int horizon,
                             std::ostream& quality_csv, Variants&& variants)
{
    const std::vector<int> successes = convergence_successes<NX, NU>(dynamics, config, variants);
    std::size_t index = 0;
    variants(
        [&](auto, char const* family, char const* algorithm, bool, bool counted)
        {
            const std::size_t slot = index++;
            if(!counted)
                return;
            write_quality_csv_row(
                quality_csv, system_name, family, algorithm, "convergence", static_cast<int>(NX), horizon,
                quality_metrics{.objective = static_cast<double>(successes[slot]) / convergence_trials,
                                .max_constraint_violation = 0.0,
                                .gradient_norm = 0.0,
                                .success = true,
                                .iterations = convergence_trials,
                                .solve_time_ms = 0.0});
        });
}

}

#endif
