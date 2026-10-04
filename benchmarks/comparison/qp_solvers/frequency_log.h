#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_FREQUENCY_LOG_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_FREQUENCY_LOG_H

#include <string>
#include <fstream>
#include <string_view>

namespace ctrlpp::bench::qp_solvers
{

/// Appends the current CPU's scaling_cur_freq (kHz) with a label, one line per sample, where the
/// kernel exposes it; the run script reduces the samples to the manifest's median, minimum and maximum.
class frequency_log
{
public:
    explicit frequency_log(const std::string& path);

    void sample(std::string_view label);

private:
    std::ofstream out_;
};

}

#endif
