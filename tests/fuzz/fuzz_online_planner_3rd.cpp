#include "ctrlpp/trajectory/online_planner_3rd.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // 8 doubles: pos, target, v_max, a_max, j_max, dt, t_start, retarget = 64 bytes
    if(size < 64)
        return 0;

    double buf[8];
    std::memcpy(buf, data, 64);

    for(int i = 0; i < 8; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    double pos = buf[0];
    double target = buf[1];
    double v_max = std::abs(buf[2]);
    double a_max = std::abs(buf[3]);
    double j_max = std::abs(buf[4]);
    double dt = std::abs(buf[5]);
    double t_start = buf[6];
    double retarget = buf[7];

    // Clamp limits to positive and reasonable
    v_max = std::clamp(v_max, 1e-3, 1e6);
    a_max = std::clamp(a_max, 1e-3, 1e6);
    j_max = std::clamp(j_max, 1e-3, 1e6);
    dt = std::clamp(dt, 1e-6, 1.0);

    ctrlpp::online_planner_3rd<double> planner({
        .v_max = v_max, .a_max = a_max, .j_max = j_max});
    planner.reset(pos);
    planner.update(target);

    // Step 10 times, retarget midway
    double t = t_start;
    for(int step = 0; step < 10; ++step)
    {
        t += dt;
        if(step == 5)
            planner.update(retarget);

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
