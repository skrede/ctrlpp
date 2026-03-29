#include "ctrlpp/trajectory/online_planner_2nd.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // 7 doubles: pos, vel, target, v_max, a_max, dt, t_start = 56 bytes
    if(size < 56)
        return 0;

    double buf[7];
    std::memcpy(buf, data, 56);

    for(int i = 0; i < 7; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    double pos = buf[0];
    double vel = buf[1];
    double target = buf[2];
    double v_max = std::abs(buf[3]);
    double a_max = std::abs(buf[4]);
    double dt = std::abs(buf[5]);
    double t_start = buf[6];

    // Clamp limits to positive and reasonable
    v_max = std::clamp(v_max, 1e-3, 1e6);
    a_max = std::clamp(a_max, 1e-3, 1e6);
    dt = std::clamp(dt, 1e-6, 1.0);

    ctrlpp::online_planner_2nd<double> planner({.v_max = v_max, .a_max = a_max});
    planner.reset(pos);

    // Set initial velocity by sampling once, then update target
    planner.update(target);

    // Step 10 times
    double t = t_start;
    for(int step = 0; step < 10; ++step)
    {
        t += dt;
        auto pt = planner.sample(t);

        if(!std::isfinite(pt.position(0)))
            __builtin_trap();
        if(!std::isfinite(pt.velocity(0)))
            __builtin_trap();
        if(!std::isfinite(pt.acceleration(0)))
            __builtin_trap();
    }

    return 0;
}
