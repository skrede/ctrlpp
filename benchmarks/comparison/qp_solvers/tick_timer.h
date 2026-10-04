#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_TICK_TIMER_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_TICK_TIMER_H

#include "arm.h"
#include "alloc_sensor.h"
#include "instruction_counter.h"

#include <chrono>
#include <cstdint>

namespace ctrlpp::bench::qp_solvers
{

using tick_clock = std::chrono::steady_clock;

struct tick_sample
{
    std::int64_t update_ns;
    std::int64_t solve_ns;
    std::int64_t tick_ns;
    std::int64_t instr_update;
    std::int64_t instr_solve;
    std::int64_t instr_tick;
    std::int64_t allocs;
};

inline auto elapsed_ns(tick_clock::time_point from, tick_clock::time_point to) -> std::int64_t
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(to - from).count();
}

inline auto counted(std::int64_t from, std::int64_t to) -> std::int64_t
{
    return from == instructions_unavailable || to == instructions_unavailable ? instructions_unavailable : to - from;
}

/// tick_ns is the one interval from the start of the update to the end of the solve. The allocation
/// sensor is armed around exactly that interval; the instruction count between the phases is read
/// inside it, by one rdpmc.
inline auto timed_tick(qp_arm& arm, const qp_tick& tick, const instruction_counter& counter) -> tick_sample
{
    arm_alloc_sensor();
    const std::int64_t i0 = counter.read();
    const tick_clock::time_point t0 = tick_clock::now();
    arm.update(tick);
    const tick_clock::time_point t1 = tick_clock::now();
    const std::int64_t i1 = counter.read();
    arm.solve();
    const tick_clock::time_point t2 = tick_clock::now();
    const std::int64_t i2 = counter.read();
    const std::int64_t allocs = disarm_alloc_sensor();
    return {elapsed_ns(t0, t1), elapsed_ns(t1, t2), elapsed_ns(t0, t2), counted(i0, i1), counted(i1, i2), counted(i0, i2), allocs};
}

}

#endif
