#include "ctrlpp/trajectory/trapezoidal_trajectory.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // 5 doubles: q0, q1, v_max, a_max, eval_time = 40 bytes
    if(size < 40)
        return 0;

    double buf[5];
    std::memcpy(buf, data, 40);

    for(int i = 0; i < 5; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    double q0 = buf[0];
    double q1 = buf[1];
    double v_max = std::abs(buf[2]);
    double a_max = std::abs(buf[3]);
    double eval_t = buf[4];

    // Clamp limits to positive
    v_max = std::clamp(v_max, 1e-6, 1e6);
    a_max = std::clamp(a_max, 1e-6, 1e6);

    ctrlpp::trapezoidal_trajectory<double> traj({
        .q0 = q0, .q1 = q1, .v_max = v_max, .a_max = a_max});

    auto pt = traj.evaluate(eval_t);

    if(!std::isfinite(pt.position(0)))
            return 0;
    if(!std::isfinite(pt.velocity(0)))
            return 0;
    if(!std::isfinite(pt.acceleration(0)))
            return 0;

    return 0;
}
