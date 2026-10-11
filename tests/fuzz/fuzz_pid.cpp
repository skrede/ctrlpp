#include "ctrlpp/control/pid.h"
#include "ctrlpp/control/pid_policies.h"

#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iterator>
#include <algorithm>

namespace
{

using PidType = ctrlpp::pid<double, 1, ctrlpp::anti_windup<ctrlpp::clamping>>;
using LimitedPid = ctrlpp::pid<double, 1, ctrlpp::rate_limit, ctrlpp::anti_windup<ctrlpp::clamping>>;
using step_result = ctrlpp::expected<PidType::vector_t, ctrlpp::pid_step_error>;

constexpr std::int32_t steps = 20;
constexpr double gain_bound = 100.0;
constexpr double signal_bound = 100.0;
constexpr double output_bound = 1e3;

// splitmix64, Steele, Lea and Flood, "Fast splittable pseudorandom number generators", OOPSLA 2014.
auto splitmix64(std::uint64_t& state) -> std::uint64_t
{
    state += 0x9e3779b97f4a7c15ULL;
    std::uint64_t z = state;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

auto stream_setpoint(std::uint64_t& state) -> double
{
    const double unit = static_cast<double>(splitmix64(state) >> 11) * 0x1.0p-53;
    return unit * 2.0 * signal_bound - signal_bound;
}

// Twice a bound on |u_k - u_{k-1}|: the proportional term kp*(sp - meas), the integrator after at most
// `steps` increments ki*e*dt (clamping only removes increments), and the previous command inside the output
// limits. The derivative is zero because the measurement is constant. The factor two keeps the rounding in
// forming the terms and rate_max*dt from closing the gap. Not finite or not positive outside the dt range the
// differential is run on.
auto idle_rate_max(double dt) -> double
{
    const double proportional = gain_bound * 2.0 * signal_bound;
    const double integral = static_cast<double>(steps) * gain_bound * 2.0 * signal_bound * dt;
    return 2.0 * (proportional + integral + output_bound) / dt;
}

// The decoded gains and limits are finite, so the only refusal the factory may
// give is the one for a lower output limit above the upper one.
auto checked_create(const PidType::config_type& cfg) -> ctrlpp::expected<PidType, ctrlpp::pid_config_error>
{
    auto created = PidType::create(cfg);
    const bool expect_refusal = cfg.output_min[0] > cfg.output_max[0];
    if(created.has_value() == expect_refusal)
        abort();
    if(!created && created.error() != ctrlpp::pid_config_error::invalid_output_limits)
        abort();
    return created;
}

auto idle_twin(const PidType::config_type& cfg, double rate_max) -> LimitedPid
{
    LimitedPid::config_type twin_cfg;
    twin_cfg.kp = cfg.kp;
    twin_cfg.ki = cfg.ki;
    twin_cfg.kd = cfg.kd;
    twin_cfg.output_min = cfg.output_min;
    twin_cfg.output_max = cfg.output_max;
    twin_cfg.policy<ctrlpp::rate_limit>().rate_max = {rate_max};
    auto created = LimitedPid::create(twin_cfg);
    if(!created)
        abort();
    return *std::move(created);
}

// The measurement is constant and there is no feed-forward, so the derivative is zero on every cycle and the
// stored integrator must reproduce the emitted command with the proportional term alone.
auto reproduced_command(const PidType& controller, double sp, double meas) -> double
{
    const auto& cfg = controller.params();
    const double unconstrained = cfg.kp[0] * (cfg.b[0] * sp - meas) + controller.integral()[0];
    return std::min(std::max(unconstrained, cfg.output_min[0]), cfg.output_max[0]);
}

// compute rejects a non-finite result in the same cycle, so a command it
// returns as a success is finite by contract.
auto checked_step(PidType& controller, const PidType::vector_t& sp, const PidType::vector_t& meas, double dt) -> step_result
{
    auto step = controller.compute(sp, meas, dt);
    const bool timestep_refused = !step && step.error() == ctrlpp::pid_step_error::invalid_timestep;
    if(timestep_refused != (dt <= 0.0))
        abort();
    if(step && !std::isfinite((*step)(0)))
        abort();
    if(step && (*step)(0) != reproduced_command(controller, sp(0), meas(0)))
        abort();
    return step;
}

auto bits(double value) -> std::uint64_t { return std::bit_cast<std::uint64_t>(value); }

auto checked_twin_step(PidType& plain, LimitedPid& twin, const PidType::vector_t& sp, const PidType::vector_t& meas, double dt) -> bool
{
    const auto plain_step = checked_step(plain, sp, meas, dt);
    const auto twin_step = twin.compute(sp, meas, dt);
    if(plain_step.has_value() != twin_step.has_value())
        abort();
    if(plain_step && bits((*plain_step)(0)) != bits((*twin_step)(0)))
        abort();
    if(!plain_step && plain_step.error() != twin_step.error())
        abort();
    if(bits(plain.integral()(0)) != bits(twin.integral()(0)))
        abort();
    return plain_step.has_value();
}

auto decode_config(const double (&buf)[8]) -> PidType::config_type
{
    PidType::config_type cfg;
    cfg.kp[0] = std::clamp(buf[0], -gain_bound, gain_bound);
    cfg.ki[0] = std::clamp(buf[1], -gain_bound, gain_bound);
    cfg.kd[0] = std::clamp(buf[2], -gain_bound, gain_bound);
    cfg.output_min[0] = std::clamp(buf[3], -output_bound, output_bound);
    cfg.output_max[0] = std::clamp(buf[4], -output_bound, output_bound);
    return cfg;
}

// The differential runs only where idle_rate_max is representable; elsewhere the controller runs alone.
void run_cycles(PidType& plain, const PidType::config_type& cfg, double setpoint, std::uint64_t stream, const PidType::vector_t& meas, double dt)
{
    const double rate_max = idle_rate_max(dt);
    if(!std::isfinite(rate_max) || rate_max <= 0.0)
    {
        for(std::int32_t step = 0; step < steps && checked_step(plain, PidType::vector_t::Constant(setpoint), meas, dt); ++step)
            setpoint = stream_setpoint(stream);
        return;
    }

    auto twin = idle_twin(cfg, rate_max);
    for(std::int32_t step = 0; step < steps && checked_twin_step(plain, twin, PidType::vector_t::Constant(setpoint), meas, dt); ++step)
        setpoint = stream_setpoint(stream);
}

}

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // 8 doubles: kp, ki, kd, output_min, output_max, dt, first setpoint, measurement = 64 bytes,
    // then an 8-byte seed for the setpoints of the later cycles.
    if(size < 72)
        return 0;

    double buf[8];
    std::memcpy(buf, data, 64);
    std::uint64_t stream;
    std::memcpy(&stream, data + 64, 8);

    if(!std::all_of(std::begin(buf), std::end(buf), [](double value) { return std::isfinite(value); }))
        return 0;

    const auto cfg = decode_config(buf);
    const double dt = buf[5];
    const PidType::vector_t meas = PidType::vector_t::Constant(std::clamp(buf[7], -signal_bound, signal_bound));
    auto created = checked_create(cfg);
    if(!created)
        return 0;

    run_cycles(*created, cfg, std::clamp(buf[6], -signal_bound, signal_bound), stream, meas, dt);
    return 0;
}
