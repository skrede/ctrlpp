#include "ctrlpp/trajectory/smoothing_spline.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // 12 doubles: 5 x + 5 y + mu + eval_point = 96 bytes
    if(size < 96)
        return 0;

    double buf[12];
    std::memcpy(buf, data, 96);

    for(int i = 0; i < 12; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    std::vector<double> x(5);
    std::vector<double> y(5);
    for(int i = 0; i < 5; ++i)
    {
        x[i] = std::clamp(buf[i], -1e6, 1e6);
        y[i] = std::clamp(buf[i + 5], -1e6, 1e6);
    }
    double mu = buf[10];
    double eval_t = buf[11];

    // Sort x values for strictly ascending with reasonable minimum span
    std::sort(x.begin(), x.end());
    for(std::size_t i = 1; i < x.size(); ++i)
    {
        if(x[i] <= x[i - 1] + 1e-3)
            x[i] = x[i - 1] + 1e-3;
    }

    // Clamp mu to valid range -- avoid extreme regularization
    mu = std::clamp(mu, 1e-3, 1.0);

    ctrlpp::smoothing_spline<double> spline({.times = x, .positions = y, .mu = mu});

    auto pt = spline.evaluate(eval_t);

    if(!std::isfinite(pt.position(0)))
        __builtin_trap();
    if(!std::isfinite(pt.velocity(0)))
        __builtin_trap();
    if(!std::isfinite(pt.acceleration(0)))
        __builtin_trap();

    return 0;
}
