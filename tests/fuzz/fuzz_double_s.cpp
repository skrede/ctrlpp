#include "ctrlpp/trajectory/double_s_trajectory.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // 6 doubles: q0, q1, v_max, a_max, j_max, eval_time = 48 bytes
    if(size < 48)
        return 0;

    double buf[6];
    std::memcpy(buf, data, 48);

    for(int i = 0; i < 6; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    double q0 = buf[0];
    double q1 = buf[1];
    double v_max = std::abs(buf[2]);
    double a_max = std::abs(buf[3]);
    double j_max = std::abs(buf[4]);
    double eval_t = buf[5];

    // Clamp limits to positive
    v_max = std::clamp(v_max, 1e-6, 1e6);
    a_max = std::clamp(a_max, 1e-6, 1e6);
    j_max = std::clamp(j_max, 1e-6, 1e6);

    ctrlpp::double_s_trajectory<double> traj({
        .q0 = q0, .q1 = q1, .v_max = v_max, .a_max = a_max, .j_max = j_max});

    auto pt = traj.evaluate(eval_t);

    if(!std::isfinite(pt.position(0)))
        __builtin_trap();
    if(!std::isfinite(pt.velocity(0)))
        __builtin_trap();
    if(!std::isfinite(pt.acceleration(0)))
        __builtin_trap();

    return 0;
}
