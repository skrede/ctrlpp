#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_RECORD_WRITER_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_RECORD_WRITER_H

#include <string>
#include <cstdint>
#include <fstream>
#include <string_view>

namespace ctrlpp::bench::qp_solvers
{

struct tick_record
{
    std::string_view cohort;
    std::string_view cell;
    std::string_view variant;
    std::int32_t nx;
    std::int32_t nu;
    std::int32_t N;
    std::int64_t n;
    std::int64_t m;
    std::string_view regime;
    std::string_view tier;
    std::string_view arm;
    std::int32_t polish;
    std::int32_t rep;
    std::int32_t seed;
    std::int32_t step;
    std::int64_t setup_ns;
    std::int64_t update_ns;
    std::int64_t solve_ns;
    std::int64_t tick_ns;
    std::int64_t condense_ns;
    std::int64_t instr_setup;
    std::int64_t instr_update;
    std::int64_t instr_solve;
    std::int64_t instr_tick;
    std::int32_t iterations;
    double kkt_stat;
    double kkt_prim;
    double kkt_comp;
    double kkt_max;
    std::string_view native_status;
    std::int32_t pass;
    std::string_view fail_class;
    std::int64_t allocs_tick;
    std::string_view dual_kind;
    std::int32_t fallback;
    double binding_frac;
    double churn_frac;
};

/// Integers in plain decimal, floats as %.17g with NaN written as nan.
class record_writer
{
public:
    record_writer(const std::string& path, std::string_view run_id);

    void write(const tick_record& row);

    auto good() const -> bool
    {
        return out_.good();
    }

private:
    std::ofstream out_;
    std::string run_id_;
};

}

#endif
