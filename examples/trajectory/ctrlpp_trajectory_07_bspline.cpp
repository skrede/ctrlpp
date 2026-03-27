// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./ctrlpp_trajectory_07_bspline' using 1:2 with lines title 'position', '' using 1:3 with lines title 'velocity'"
// Redirect: ./ctrlpp_trajectory_07_bspline > output.csv

#include <ctrlpp/trajectory/bspline_trajectory.h>

#include <iomanip>
#include <iostream>

int main()
{
    // CNC tool-tip path through 6 control points using cubic B-spline (Degree=3).
    // Empty knot vector triggers auto-generated uniform clamped knots.
    ctrlpp::bspline_trajectory<double, 3>::config cfg{
        .control_points = {0.0, 2.0, 5.0, 4.0, 7.0, 10.0},
    };

    ctrlpp::bspline_trajectory<double, 3> bspline(cfg);

    auto const total = bspline.duration();
    constexpr double dt = 0.01;

    std::cout << "time,position,velocity\n";

    for (double t = 0.0; t <= total; t += dt)
    {
        auto const pt = bspline.evaluate(t);
        std::cout << std::fixed << std::setprecision(4) << t << "," << pt.position[0] << "," << pt.velocity[0] << "\n";
    }
}
