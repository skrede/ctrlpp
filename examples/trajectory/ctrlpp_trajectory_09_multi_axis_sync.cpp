// Usage: ./ctrlpp_trajectory_09_multi_axis_sync | gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '-' using 1:2 with lines title 'x_pos', '' using 1:4 with lines title 'y_pos', '' using 1:6 with lines title 'z_pos'"
// Redirect: ./ctrlpp_trajectory_09_multi_axis_sync > output.csv

#include <ctrlpp/trajectory/synchronize.h>
#include <ctrlpp/trajectory/trapezoidal_trajectory.h>

#include <algorithm>
#include <iomanip>
#include <iostream>

int main()
{
    // 3-axis robot move: X 0->100mm, Y 0->50mm, Z 0->200mm.
    // All axes share v_max=80mm/s, a_max=400mm/s^2.
    // Synchronize so all axes finish simultaneously.
    ctrlpp::trapezoidal_trajectory<double> ax_x({.q0 = 0.0, .q1 = 100.0, .v_max = 80.0, .a_max = 400.0});
    ctrlpp::trapezoidal_trajectory<double> ax_y({.q0 = 0.0, .q1 = 50.0, .v_max = 80.0, .a_max = 400.0});
    ctrlpp::trapezoidal_trajectory<double> ax_z({.q0 = 0.0, .q1 = 200.0, .v_max = 80.0, .a_max = 400.0});

    ctrlpp::synchronize(ax_x, ax_y, ax_z);

    std::cerr << "Synchronized duration: " << ax_x.duration() << " s\n";

    auto const total = ax_x.duration();
    constexpr double dt = 0.001;

    std::cout << "time,x_pos,x_vel,y_pos,y_vel,z_pos,z_vel\n";

    for (double t = 0.0; t <= total; t += dt)
    {
        auto const px = ax_x.evaluate(t);
        auto const py = ax_y.evaluate(t);
        auto const pz = ax_z.evaluate(t);

        std::cout << std::fixed << std::setprecision(4) << t << "," << px.position[0] << "," << px.velocity[0] << "," << py.position[0] << "," << py.velocity[0] << "," << pz.position[0] << "," << pz.velocity[0] << "\n";
    }
}
