#include "ctrlpp/control/pid.h"
#include "ctrlpp/control/pid_policies.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <algorithm>

namespace
{

using PidType = ctrlpp::pid<double, 1, ctrlpp::anti_windup<ctrlpp::clamping>>;

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

// compute rejects a non-finite result in the same cycle, so a command it
// returns as a success is finite by contract.
auto checked_step(PidType& controller, const PidType::vector_t& sp, const PidType::vector_t& meas, double dt) -> bool
{
    auto step = controller.compute(sp, meas, dt);
    const bool timestep_refused = !step && step.error() == ctrlpp::pid_step_error::invalid_timestep;
    if(timestep_refused != (dt <= 0.0))
        abort();
    if(step && !std::isfinite((*step)(0)))
        abort();
    return step.has_value();
}

}

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // 8 doubles: kp, ki, kd, output_min, output_max, dt, setpoint, measurement = 64 bytes
    if(size < 64)
        return 0;

    double buf[8];
    std::memcpy(buf, data, 64);

    for(std::int32_t i = 0; i < 8; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    PidType::config_type cfg;
    cfg.kp[0] = std::clamp(buf[0], -100.0, 100.0);
    cfg.ki[0] = std::clamp(buf[1], -100.0, 100.0);
    cfg.kd[0] = std::clamp(buf[2], -100.0, 100.0);
    cfg.output_min[0] = std::clamp(buf[3], -1e3, 1e3);
    cfg.output_max[0] = std::clamp(buf[4], -1e3, 1e3);
    const double dt = buf[5];
    const PidType::vector_t sp = PidType::vector_t::Constant(std::clamp(buf[6], -100.0, 100.0));
    const PidType::vector_t meas = PidType::vector_t::Constant(std::clamp(buf[7], -100.0, 100.0));

    auto created = checked_create(cfg);
    if(!created)
        return 0;

    for(std::int32_t step = 0; step < 20; ++step)
    {
        if(!checked_step(*created, sp, meas, dt))
            return 0;
    }

    return 0;
}
