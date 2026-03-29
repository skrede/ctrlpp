#include "ctrlpp/trajectory/bspline_trajectory.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size)
{
    // Cubic B-spline (degree 3) with 5 control points
    // Knot vector: 5 + 3 + 1 = 9 entries
    // Total: 5 control points + 9 knots + 1 eval_t = 15 doubles = 120 bytes
    if(size < 120)
        return 0;

    double buf[15];
    std::memcpy(buf, data, 120);

    for(int i = 0; i < 15; ++i)
    {
        if(!std::isfinite(buf[i]))
            return 0;
    }

    std::vector<double> control_points(buf, buf + 5);
    std::vector<double> knots(buf + 5, buf + 14);
    double eval_t = buf[14];

    // Sort knot vector to ensure non-decreasing
    std::sort(knots.begin(), knots.end());

    // Ensure knot vector has proper clamped structure for degree 3, 5 control points
    // Required: first 4 equal, last 4 equal, interior monotone
    // Just sort and let the constructor validate
    try
    {
        ctrlpp::bspline_trajectory<double, 3> bspline({
            .control_points = control_points,
            .knot_vector = knots,
        });

        auto pt = bspline.evaluate(eval_t);

        if(!std::isfinite(pt.position(0)))
            __builtin_trap();
        if(!std::isfinite(pt.velocity(0)))
            __builtin_trap();
        if(!std::isfinite(pt.acceleration(0)))
            __builtin_trap();
    }
    catch(const std::invalid_argument&)
    {
        // Invalid B-spline config is expected for some fuzz inputs
    }

    return 0;
}
