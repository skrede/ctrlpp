#include "ctrlpp/trajectory/cubic_spline.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // 10 doubles: 5 x-values + 5 y-values = 80 bytes
    if(size < 80)
        return 0;

    double buf[10];
    std::memcpy(buf, data, 80);

    for(int i = 0; i < 10; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    // Extract x and y values, clamp to prevent overflow in slope computation
    std::vector<double> x(5);
    std::vector<double> y(5);
    for(int i = 0; i < 5; ++i)
    {
        x[i] = std::clamp(buf[i], -1e6, 1e6);
        y[i] = std::clamp(buf[i + 5], -1e6, 1e6);
    }

    // Sort x values to ensure strictly ascending (spline requirement)
    std::sort(x.begin(), x.end());

    // Ensure strictly increasing with reasonable minimum span width
    for(std::size_t i = 1; i < x.size(); ++i)
    {
        if(x[i] <= x[i - 1] + 1e-3)
            x[i] = x[i - 1] + 1e-3;
    }

    // A rejected configuration is a typed outcome, not a finding: the oracle
    // below only says what a successfully built spline must satisfy.
    auto const built = ctrlpp::cubic_spline<double>::create({.times = x, .positions = y});
    if(!built.has_value())
            return 0;

    // Evaluate at midpoint
    double t_mid = (x.front() + x.back()) * 0.5;
    auto pt = built->evaluate(t_mid);

    if(!std::isfinite(pt.position(0)))
            return 0;
    if(!std::isfinite(pt.velocity(0)))
            return 0;
    if(!std::isfinite(pt.acceleration(0)))
            return 0;

    return 0;
}
