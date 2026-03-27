// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./ctrlpp_trajectory_06_cubic_spline' using 1:2 with lines title 'position', '' using 1:3 with lines title 'velocity'"
// Redirect: ./ctrlpp_trajectory_06_cubic_spline > output.csv

#include <ctrlpp/trajectory/cubic_spline.h>

#include <iomanip>
#include <iostream>

int main()
{
    // Robot pick-and-place: joint passes through 5 waypoints at specified times.
    // Natural boundary conditions (zero acceleration at endpoints).
    ctrlpp::cubic_spline<double>::config cfg{
        .times = {0.0, 1.0, 2.5, 3.5, 5.0},
        .positions = {0.0, 1.2, 0.8, 2.0, 1.5},
        .bc = ctrlpp::boundary_condition::natural,
    };

    ctrlpp::cubic_spline<double> spline(cfg);

    constexpr double dt = 0.01;

    std::cout << "time,position,velocity,acceleration\n";

    for (double t = cfg.times.front(); t <= cfg.times.back(); t += dt)
    {
        auto const pt = spline.evaluate(t);
        std::cout << std::fixed << std::setprecision(4) << t << "," << pt.position[0] << "," << pt.velocity[0] << "," << pt.acceleration[0] << "\n";
    }
}
